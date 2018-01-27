
#include "HartreeFock.hpp"

int main() {
    
    HartreeFock HF(0.01, 0.01);
    HF.print_state();
    HF.Iterate();
    
}
