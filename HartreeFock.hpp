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

    double EMP2;

    Matrix S;
    Matrix T;
    Matrix V;
    Matrix Hcore;

    Vector E;

    Matrix SOM;
    Matrix F0;
    Matrix FMO;
    Matrix C0;
    Matrix e0;
    Matrix D0;
    Matrix prev_D0;

    TEIMatrix TEI;
    TEIMatrix TEI_MO;

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
    void MOBasisFock();
    bool EConverg();
    bool DensConverg();
    void SaveEnergy();
    void SaveDensity();

    void MP2_Correction();
    void Set_OrbitalEnergy();
    double MP2_Energy();
    void TEI_Transform_N5();
    void TEI_Transform_N8();
    void CheckEnergy();

    HartreeFock(double tol_e, double tol_dens);

};

#endif /* HartreeFock_hpp */
