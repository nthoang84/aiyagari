#ifndef AIYAGARI_H
#define AIYAGARI_H

#include <vector>

using namespace std;

class Aiyagari {
public:
    Aiyagari();

    void discretizeLabor(double mTauchen = 3.0);

    void computeLaborInvDist(double eps = EPS, int maxIter = MAX_ITER);

    void print() const;

private:
    static constexpr double EPS = 1e-6;
    static constexpr int MAX_ITER = (int) 1e4;

    int assetGridSize;
    int laborGridSize;
    double rho;
    double sigma;

    vector<double> labor;
    vector<vector<double>> transition;
    vector<double> laborInvDist;
};

#endif