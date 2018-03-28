
#include "HartreeFock.hpp"
using namespace std;

int main(int argc, char **argv) {

    string testCase = "STO3G_Water"; //"DZ_Water"; //

    HartreeFock HF1(testCase, 1e-6, 1e-6);
    HF1.Iterate();
    HF1.MP2_Correction();

    HartreeFock HF2(testCase, 1e-6, 1e-6);
    HF2.DIISIterate();
    HF2.MP2_Correction();

    return 0;

}
