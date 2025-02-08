#include "aiyagari.h"

#include <cmath>

#include <iostream>

Aiyagari::Aiyagari() :
    assetGridSize(300),             
    laborGridSize(7),               // Aiyagari (1994, p. 675)
    rho(0.4),                       // [0.2, 0.4] as in Aiyagari (1994, p. 675)
    sigma(0.9)                      // {0, 0.3, 0.6, 0.9} as in Aiyagari (1994, p. 675)
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
        return std::erfc(-x / std::sqrt(2)) / 2;
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

    for (int i = 0; i < laborGridSize; i++) {
        double sum = 0;
        for (int j = 0; j < laborGridSize; j++) {
            cout << transition[i][j] << ' ';
            sum += transition[i][j];
        }
        cout << "\tSUM:" << sum;
        cout << '\n';
    }
}

void Aiyagari::print() const {
    cout << assetGridSize << ' ' << laborGridSize << endl;
}