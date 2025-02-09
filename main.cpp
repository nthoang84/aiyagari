#include "aiyagari.h"

int main() {
    Aiyagari model;
    model.discretizeLabor();
    model.computeLaborInvDist();
    model.computeAssetGrid();
    model.solveEquilibrium();
    model.print();
    return 0;
}