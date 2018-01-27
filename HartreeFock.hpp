//
//  HartreeFock.hpp
//  HartreeFock
//
//  Created by Isabel Craig on 1/26/18.
//  Copyright © 2018 Isabel Craig. All rights reserved.
//

#ifndef HartreeFock_hpp
#define HartreeFock_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

#include "Read.hpp"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

void Diagonlize(Matrix *M, Matrix *evals, Matrix *evecs);

class HartreeFock {

public:

    double enuc;
    double tol_dens;
    double tol_e;
    
    Matrix S;
    Matrix T;
    Matrix V;
    Matrix Hcore;

    Matrix SOM;
    Matrix F0;
    Matrix C0;
    Matrix e0;
    Matrix D0;
    Matrix prev_D0;
    MullikenMatrix TEI;
    
    double eelec;
    double etot;
    double prev_etot;
    double delE;
    double rmsD;

    void print_state();
    void Set_DensityMatrix();
    void Set_InitialFock();
    void Set_Fock();
    void SymmetricOrth();
    void Iterate();
    void Set_Energy();
    bool EConverg();
    bool DensConverg();

    HartreeFock(double tol_e, double tol_dens);
    
};

#endif /* HartreeFock_hpp */
