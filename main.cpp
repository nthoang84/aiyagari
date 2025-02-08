#include "aiyagari.h"

int main() {
    Aiyagari model;
    model.discretizeLabor();
    model.computeLaborInvDist();
    return 0;
}