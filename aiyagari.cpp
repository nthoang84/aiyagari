#include "aiyagari.h"

#include <cmath>

#include <iostream>

Aiyagari::Aiyagari() 
    : assetGridSize(300),
      laborGridSize(7),               // Aiyagari (1994, p. 675)
      alpha(0.36),                    // Capital share of income
      beta(0.96),                     // Discount factor
      gamma(3.0),                     // Coefficient of relative risk aversion
      rho(0.4),                       // [0.2, 0.4] as in Aiyagari (1994, p. 675)
      sigma(0.9),                     // {0, 0.3, 0.6, 0.9} as in Aiyagari (1994, p. 675)
      depreciationRate(0.05),
      borrowingLimit(0.0),
      assetMax(500)
{
    totalGridSize = assetGridSize * laborGridSize;
    assetMin = -borrowingLimit;
}

void Aiyagari::discretizeLabor(double mTauchen) {
    // Labor process is discretized using Tauchen (1986) method.
    double varShock = (sigma * sigma) / (1 - rho * rho);
    double stdShock = sqrt(varShock);
    labor.resize(laborGridSize);
    labor[laborGridSize - 1] = mTauchen * stdShock;
    labor[0] = -labor[laborGridSize - 1];
    double stepSize = (labor.back() - labor.front()) / (laborGridSize - 1);
    for (int i = 1; i < laborGridSize - 1; i++) {
        labor[i] = labor[i - 1] + stepSize;
    }

    // Method to calculate standard normal CDF
    auto phi = [&](double x) {
        return erfc(-x / sqrt(2)) / 2;
    };

    // Compute Markov transition probability matrix
    transition.resize(laborGridSize);
    for (int i = 0; i < laborGridSize; i++) {
        transition[i].resize(laborGridSize);
        for (int j = 0; j < laborGridSize; j++) {
            double left = (labor[j] - rho * labor[i] - stepSize / 2) / sigma;
            double right = (labor[j] - rho * labor[i] + stepSize / 2) / sigma;
            if (j == 0) {
                transition[i][j] = phi(right);
            } else if (j == laborGridSize - 1) {
                transition[i][j] = 1 - phi(left);
            } else {
                transition[i][j] = phi(right) - phi(left);
            } 
        }
    }

    // Exponentiate to get labor from log(labor)
    aggregateLabor = 0;
    for (int i = 0; i < laborGridSize; i++) {
        labor[i] = exp(labor[i]);
        aggregateLabor += labor[i];
    }
}

void Aiyagari::computeLaborInvDist(double eps, int maxIter) {
    laborInvDist.resize(laborGridSize);
    laborInvDist[0] = 1.0;
    int iter = 0;
    while (iter < maxIter) {
        vector<double> dist(laborGridSize);
        for (int j = 0; j < laborGridSize; j++) {
            for (int i = 0; i < laborGridSize; i++) {
                dist[j] += laborInvDist[i] * transition[i][j];
            }
        }
        double diff = 0;
        for (int i = 0; i < laborGridSize; i++) {
            diff = max(diff, fabs(dist[i] - laborInvDist[i]));
            laborInvDist[i] = dist[i];
        }
        if (diff < eps) {
            break;
        }
        iter++;
    }
    double sum = 0;
    for (int i = 0; i < laborGridSize; i++) {
        sum += laborInvDist[i];
    }
    if (fabs(sum - 1.0) > EPS) {
        for (int i = 0; i < laborGridSize; i++) {
            laborInvDist[i] /= sum;
        }
    }
}

void Aiyagari::computeAssetGrid(double growthRate) {
    asset.resize(assetGridSize);
    if (fabs(growthRate) < EPS) {    
        double stepSize = (assetMax - assetMin) / (assetGridSize - 1);
        for (int i = 0; i < assetGridSize; i++) {
            asset[i] = (i == 0 ? assetMin : asset[i - 1] + stepSize);
        }
        return;
    }
    for (int i = 0; i < assetGridSize; i++) {
        asset[i] = assetMin + (assetMax - assetMin) * 
                   ((pow(1 + growthRate, i) - 1) / (pow(1 + growthRate, assetGridSize - 1) - 1));
    }
}

void Aiyagari::computePolicy(double interestRate, double eps, int maxIter) {
    double wageRate = computeWageFromInterestRate(interestRate);
    assetPolicy.resize(totalGridSize);
    consumptionPolicy.resize(totalGridSize);
    endogenousAsset.resize(totalGridSize);
    MU.resize(totalGridSize);
    expectedMU.resize(totalGridSize);

    // Given the policy function for tomorrow's asset, calculate all relevant quantities 
    // in the Euler equation and update the endogenous asset grid.
    auto updateGrid = [&]() {
        for (int j = 0; j < laborGridSize; j++) {
            for (int i = 0; i < assetGridSize; i++) {
                consumptionPolicy[id(i, j)] = (1 + interestRate) * asset[i] + 
                                              wageRate * labor[j] -
                                              assetPolicy[id(i, j)];
                MU[id(i, j)] = mu_c(consumptionPolicy[id(i, j)]);
            }
        }
        for (int j = 0; j < laborGridSize; j++) {
            for (int i = 0; i < assetGridSize; i++) {
                expectedMU[id(i, j)] = 0;
                for (int k = 0; k < laborGridSize; k++) {
                    expectedMU[id(i, j)] += beta * (1 + interestRate) * 
                                           transition[j][k] * MU[id(i, k)];
                }
                endogenousAsset[id(i, j)] = (mu_c_inverse(expectedMU[id(i, j)]) + 
                                             asset[i] - wageRate * labor[j]) / (1 + interestRate);
            }
        }
    };
    
    // Start with an initial guess where tomorrow's assets are zero for all scenarios.
    updateGrid();

    int iter = 0;
    while (iter < maxIter) {       
        for (int j = 0; j < laborGridSize; j++) {
            int current_i = 0;
            for (int i = 0; i < assetGridSize; i++) {
                while (current_i < assetGridSize && endogenousAsset[id(current_i, j)] < asset[i]) {
                    current_i++;
                }
                if (current_i == 0) {
                    assetPolicy[id(i, j)] = asset[0];
                } else if (current_i == assetGridSize) {
                    assetPolicy[id(i, j)] = asset.back();
                } else {
                    const double& left = endogenousAsset[id(current_i - 1, j)];
                    const double& right = endogenousAsset[id(current_i, j)];
                    double weight = (asset[i] - left) / (right - left);
                    assetPolicy[id(i, j)] = (1 - weight) * asset[current_i - 1] + weight * asset[current_i];
                }
            }
        }
        auto prevMU = MU;
        updateGrid();
        double diff = 0;
        for (int i = 0; i < assetGridSize; i++) {
            for (int j = 0; j < laborGridSize; j++) {
                diff = max(diff, fabs(MU[id(i, j)] - prevMU[id(i, j)]));
            }
        }
        if (diff < eps) {
            break;
        }
        iter++;
    }
}

void Aiyagari::simulate(double eps, int maxIter) {
    // Approximate the density using a histogram over a fixed grid (Young, 2010, JEDC).
    vector<double> where(totalGridSize);
    vector<double> weight(totalGridSize);
    for (int j = 0; j < laborGridSize; j++) {
        int current_i = 0;
        for (int i = 0; i < assetGridSize; i++) {
            double x = assetPolicy[id(i, j)];
            while (current_i < assetGridSize && asset[current_i] < x) {
                current_i++;
            }
            where[id(i, j)] = current_i;
            if (current_i == 0 || current_i == assetGridSize) {
                continue;
            } else {
                weight[id(i, j)] = (x - asset[current_i - 1]) / (asset[current_i] - asset[current_i - 1]);
            }
        }
    }
    vector<double> dist(totalGridSize, 1.0 / totalGridSize);
    int iter = 0;
    while (iter < maxIter) {
        vector<double> newDist(totalGridSize);
        for (int j = 0; j < laborGridSize; j++) {
            for (int i = 0; i < assetGridSize; i++) {
                for (int k = 0; k < laborGridSize; k++) {
                    double new_i = where[id(i, j)];
                    double w = weight[id(i, j)];
                    if (new_i == 0) {
                        newDist[id(0, k)] += dist[id(i, j)] * transition[j][k];
                    } else if (new_i == assetGridSize) {
                        newDist[id(assetGridSize - 1, k)] += dist[id(i, j)] * transition[j][k];
                    } else {
                        newDist[id(new_i - 1, j)] += dist[id(i, j)] * transition[j][k] * (1 - w);
                        newDist[id(new_i, j)] += dist[id(i, j)] * transition[j][k] * w;
                    }
                }
            }
        }
        double diff = 0;
        for (int i = 0; i < assetGridSize; i++) {
            for (int j = 0; j < laborGridSize; j++) {
                diff = max(diff, fabs(dist[id(i, j)] - newDist[id(i, j)]));
            }
        }
        for (int i = 0; i < totalGridSize; i++) {
            dist[i] = newDist[i];
        }
        if (diff < eps) {
            break;
        }
        iter++;
    }
    aggregateCapitalSupply = 0;
    for (int i = 0; i < assetGridSize; i++) {
        for (int j = 0; j < laborGridSize; j++) {
            aggregateCapitalSupply += asset[i] * dist[id(i, j)];
        }
    }
}

void Aiyagari::solveEquilibrium(double eps, int maxIter) {
    double interestRate = 0.04;
    int iter = 0;
    while (iter < maxIter) {
        double wageRate = computeWageFromInterestRate(interestRate);
        aggregateCapitalDemand = aggregateLabor * pow(alpha / (interestRate + depreciationRate), 1.0 / (1.0 - alpha));
        computePolicy(interestRate);
        simulate();        
        double diff = (aggregateCapitalSupply - aggregateCapitalDemand) / 
                      ((aggregateCapitalSupply + aggregateCapitalDemand) / 2);
        if (fabs(diff) < eps) {
            break;
        }
        interestRate = alpha * pow(aggregateLabor, 1 - alpha) / pow((aggregateCapitalSupply + aggregateCapitalDemand) / 2, 1 - alpha) - depreciationRate;
        iter++;
        cout << "Iteration " << iter << ": excessDemand = " << aggregateCapitalDemand - aggregateCapitalSupply << '\n';
    }
    eqmInterestRate = interestRate;
}

void Aiyagari::print() const {
    cout << ">> Equilibrium: r = " << eqmInterestRate << ", K(r) = " << aggregateCapitalDemand << '\n';
}

void Aiyagari::plot() {
    pair<double, double> interestRateBounds{0.005, 0.08};
    double interestRateStep = 0.0025;
    vector<double> interestRate;
    for (double currentRate = interestRateBounds.first; currentRate < interestRateBounds.second; ) {
        interestRate.push_back(currentRate);
        currentRate += interestRateStep;
        if (fabs(currentRate - interestRateBounds.second) < EPS) {
            interestRate.push_back(currentRate);
        }
    }
    vector<double> demand(interestRate.size());
    vector<double> supply(interestRate.size());
    for (int i = 0; i < interestRate.size(); i++) {
        double wageRate = computeWageFromInterestRate(interestRate[i]);
        demand[i] = aggregateLabor * pow(alpha / (interestRate[i] + depreciationRate), 1.0 / (1.0 - alpha));
        computePolicy(interestRate[i]);
        simulate();
        supply[i] = aggregateCapitalSupply;
    }
    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (!gnuplot) {
        std::cerr << "Error: Unable to open gnuplot." << std::endl;
        return;
    }
    fprintf(gnuplot, "set terminal pdfcairo\n");
    fprintf(gnuplot, "set output 'capital_supply_demand.pdf'\n");
    fprintf(gnuplot, "set xlabel 'Aggregate capital'\n");
    fprintf(gnuplot, "set ylabel 'Interest rate'\n");
    fprintf(gnuplot, "plot '-' with lines title 'Demand', '-' with lines title 'Supply'\n");
    for (size_t i = 0; i < interestRate.size(); ++i) {
        fprintf(gnuplot, "%lf %lf\n", demand[i], interestRate[i]);
    }
    fprintf(gnuplot, "e\n");
    for (size_t i = 0; i < interestRate.size(); ++i) {
        fprintf(gnuplot, "%lf %lf\n", supply[i], interestRate[i]);
    }
    fprintf(gnuplot, "e\n");
    fclose(gnuplot);
}

inline int Aiyagari::id(int x, int y) {
    return x * laborGridSize + y;
}

inline double Aiyagari::u(double c) {
    return pow(c, 1 - gamma) / (1 - gamma);
}

inline double Aiyagari::mu_c(double c) {
    return pow(c, -gamma);
}

inline double Aiyagari::mu_c_inverse(double u) {
    return pow(u, -1 / gamma);
}

inline double Aiyagari::computeWageFromInterestRate(double interestRate) {
    return (1 - alpha) * pow(alpha / (interestRate + depreciationRate), alpha / (1 - alpha));
}