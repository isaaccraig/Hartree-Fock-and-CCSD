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
#include "QuantumUtils.hpp"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

class HartreeFock {

  public:

    double enuc;
    double tol_dens;
    double tol_e;
    double EMP2;

    Eigen::MatrixXd *S;
    Eigen::MatrixXd *V;
    Eigen::MatrixXd *T;
    Eigen::MatrixXd *Hcore;
    Eigen::MatrixXd *E;
    Eigen::MatrixXd *SOM;
    Eigen::MatrixXd *F0;
    Eigen::MatrixXd *FMO;
    Eigen::MatrixXd *C0;
    Eigen::MatrixXd *e0;
    Eigen::MatrixXd *D0;
    Eigen::MatrixXd *prev_D0;

    Eigen::MatrixXd *TEI_MO;
    Eigen::MatrixXd *TEI_AO;

    double eelec;
    double etot;
    double prev_etot;
    double delE;
    double rmsD;

    void print_state();
    void Set_DensityMatrix();
    void Set_InitialFock();
    void Set_Fock();
    void Iterate();
    void Set_Energy();
    void Set_MOBasisFock();
    bool EConverg();
    bool DensConverg();
    void SaveEnergy();
    void SaveDensity();

    void MP2_Correction();
    void Set_OrbitalEnergy();
    void CheckEnergy();

    HartreeFock(double tol_e, double tol_dens);
    ~HartreeFock();

};

#endif /* HartreeFock_hpp */
