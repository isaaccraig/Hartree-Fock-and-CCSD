
#include "HartreeFock.hpp"
using namespace std;

int main() {

    try {
        HartreeFock HF(0.001, 0.001);
        HF.print_state();
        HF.Iterate();
        HF.MP2_Correction();
        HF.CheckEnergy();
    } catch (SanityCheckException e) {
        cout << e.what() << endl;
    }

    return 0;

}
