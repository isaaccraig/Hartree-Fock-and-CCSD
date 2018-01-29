
#include "HartreeFock.hpp"

int main() {

    HartreeFock HF(0.01, 0.01);
    HF.PrintState();
    HF.Iterate( 2 );

}
