#include "aiyagari.h"

#include <iostream>
#include <chrono>

int main() {
    auto start = chrono::high_resolution_clock::now();
    Aiyagari model;
    model.discretizeLabor();
    model.computeLaborInvDist();
    model.computeAssetGrid();
    model.solveEquilibrium();
    model.print();
    model.plot();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    cout << ">> Elapsed time: " << duration.count() << " seconds" << endl;
    return 0;
}