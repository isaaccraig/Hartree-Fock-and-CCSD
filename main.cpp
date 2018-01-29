
#include "HartreeFock.hpp"

int main() {
    
    HartreeFock HF(0.001, 0.001);
    HF.print_state();
    HF.Iterate();
    
    
    
}
