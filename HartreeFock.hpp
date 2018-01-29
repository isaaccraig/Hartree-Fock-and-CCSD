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

#include "Utils.hpp"

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
    TEIMatrix TEI;

    double eelec;
    double etot;
    double prev_etot;
    double delE;
    double rmsD;
    double mu_x;
    double mu_y;
    double mu_z;
    double tot_dip_moment;
    double q0;
    double q1;
    double q2;

    void PrintState();
    void SaveDensity();
    void SaveEnergy();
    void SetDensityMatrix();
    void SetInitialFock();
    void SetFock();
    void Iterate( int maxit );
    void SetEnergy();
    bool EConverg();
    bool DensConverg();

    HartreeFock(double tol_e, double tol_dens);

};

#endif /* HartreeFock_hpp */
