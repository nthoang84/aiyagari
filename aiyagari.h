#ifndef AIYAGARI_H
#define AIYAGARI_H

#include <vector>

using namespace std;

class Aiyagari {
public:
    Aiyagari();

    void discretizeLabor(double mTauchen = 3.0);

    void print() const;

private:
    int assetGridSize;
    int laborGridSize;
    double rho;
    double sigma;

    vector<double> labor;
    vector<vector<double>> transition;
};

#endif