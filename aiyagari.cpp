#include "aiyagari.h"

#include <cmath>

#include <iostream>

Aiyagari::Aiyagari() :
    assetGridSize(300),             
    laborGridSize(7),               // Aiyagari (1994, p. 675)
    rho(0.4),                       // [0.2, 0.4] as in Aiyagari (1994, p. 675)
    sigma(0.9),                     // {0, 0.3, 0.6, 0.9} as in Aiyagari (1994, p. 675)
    borrowingLimit(0.0),            
    assetMax(300.0)                 
{}

/**
 * Labor endowments follow a Markov chain as in Aiyagari (1994, p. 675).
 * Labor process is discretized using Tauchen (1986) method.
 */
void Aiyagari::discretizeLabor(double mTauchen) {
    // Discretize labor process
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
    assetMin = -borrowingLimit;
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

void Aiyagari::print() const {
    cout << assetGridSize << ' ' << laborGridSize << endl;
}