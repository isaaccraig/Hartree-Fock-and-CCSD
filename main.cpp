
#include "HartreeFock.hpp"
#include "CoupledCluster.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv) {

    string testCase = "STO3G_Water"; //"DZ_Water"; //

    HartreeFock HF1(testCase, 1e-6, 1e-6);
    HF1.Iterate();
    HF1.MP2_Correction();

    HartreeFock HF2(testCase, 1e-6, 1e-6);
    HF2.DIISIterate();
    HF2.MP2_Correction();

    CCSD ccsd(HF1, 1e-6);
    ccsd.Iterate();
    double mp2 = ccsd.MP2Energy();
    double corrE = ccsd.correlationEnergy();

    return 0;

}
