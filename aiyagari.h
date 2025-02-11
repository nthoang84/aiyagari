#ifndef AIYAGARI_H
#define AIYAGARI_H

#include <vector>

using namespace std;

class Aiyagari {
public:
    Aiyagari();

    void discretizeLabor(double mTauchen = 3.0);

    void computeLaborInvDist(double eps = EPS, int maxIter = MAX_ITER);

    void computeAssetGrid(bool plotDistribution = false, double growthRate = 0.025);

    void computePolicy(double interestRate, double eps = EPS, int maxIter = MAX_ITER);

    void simulate(double eps = EPS, int maxIter = MAX_ITER);

    void solveEquilibrium(double eps = EPS, int maxIter = MAX_ITER);

    void print() const;

    void plot(bool verbose = false);

    void plot(const vector<double>& data, const string& label, bool isGrid = false);

private:
    inline int id(int x, int y);
    inline double u(double c);
    inline double mu_c(double c);
    inline double mu_c_inverse(double u);
    inline double computeWageFromInterestRate(double interestRate);

    static constexpr double EPS = 1e-6;
    static constexpr int MAX_ITER = (int) 1e3;

    int assetGridSize;
    int laborGridSize;
    int totalGridSize;
    double alpha;
    double beta;
    double gamma;
    double rho;
    double sigma;
    double depreciationRate;
    double borrowingLimit;
    double assetMin;
    double assetMax;
    double aggregateCapitalSupply;
    double aggregateCapitalDemand;
    double aggregateLabor;
    double eqmInterestRate;

    vector<double> labor;
    vector<vector<double>> transition;
    vector<double> laborInvDist;
    vector<double> asset;
    vector<double> assetPolicy;
    vector<double> consumptionPolicy;
    vector<double> endogenousAsset;
    vector<double> endogenousAssetBound;
    vector<double> MU;
    vector<double> expectedMU;
};

#endif