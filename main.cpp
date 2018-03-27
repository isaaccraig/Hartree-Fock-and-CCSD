
#include "HartreeFock.hpp"
using namespace std;

int main(int argc, char **argv) {

    string testCase = "STO3G_Water"; //"DZ_Water"; //

    HartreeFock HF(testCase, 1e-3, 1e-3);
    HF.print_state();
    HF.Iterate();
    HF.MP2_Correction();
    HF.CheckEnergy();
    HF.DipoleMoment();

    return 0;

}
