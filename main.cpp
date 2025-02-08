#include "aiyagari.h"

int main() {
    Aiyagari model;
    model.discretizeLabor();
    model.computeLaborInvDist();
    model.computeAssetGrid();
    return 0;
}