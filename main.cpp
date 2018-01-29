
#include "HartreeFock.hpp"

int main() {

    HartreeFock HF(0.001, 0.001);
    HF.print_state();
    HF.Iterate();
    HF.MP2_Correction();
    HF.CheckEnergy();

    return 0;

}
