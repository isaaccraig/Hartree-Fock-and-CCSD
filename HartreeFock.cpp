
//
//  main.cpp
//  HartreeFock
//
//  Created by Isabel Craig on 1/26/18.
//  Copyright © 2018 Isabel Craig. All rights reserved.
//

#include "HartreeFock.hpp"
#include "Utils.hpp"

using namespace std;

HartreeFock::HartreeFock(double tol_dens, double tol_e){

    READIN::val("data/enuc.dat", &enuc);
    READIN::SymMatrix("data/overlap.dat",  &S);
    READIN::SymMatrix("data/kinetic.dat", &T);
    READIN::SymMatrix("data/anuc.dat", &V);
    READIN::TEI("data/eri.dat", &TEI);
    UTILS::SymmetricOrth(&S, &SOM);     // Symmetric Orthogalization Matrix

    Hcore = T + V;

    SetInitialFock();            // Build Initial Guess Fock Matrix
    SetDensityMatrix();          // Build Initial Density Matrix using occupied MOs
    SetEnergy();                 // Compute the Initial SCF Energy

}

void HartreeFock::PrintState() {

    cout << "Nuclear repulsion energy = \n" << enuc << endl;
    cout << " -------------------------------------------- " << endl;
    cout << "Overlap Integrals: \n" << S << endl;
    cout << " -------------------------------------------- " << endl;
    cout << "Kinetic-Energy Integrals: \n" << T << endl;
    cout << " -------------------------------------------- " << endl;
    cout << "Nuclear Attraction Integrals: \n" << V << endl;
    cout << " -------------------------------------------- " << endl;
    cout << "Core Hamiltonian: \n" << Hcore << endl;
    cout << " -------------------------------------------- " << endl;
    cout << "S^-1/2 Matrix: \n" << SOM << endl;
    cout << " -------------------------------------------- " << endl;
    cout << "Initial F' Matrix: \n" << F0 << endl;
    cout << " -------------------------------------------- " << endl;
    cout << "Initial C Matrix: \n" << C0 << endl;
    cout << " -------------------------------------------- " << endl;
    cout << "Initial Density Matrix: \n" << D0 << endl;
    cout << " -------------------------------------------- " << endl;
    cout << "Initial Energy: \n" << etot << endl;
    cout << " -------------------------------------------- " << endl;

}

void HartreeFock::SetInitialFock(){
    F0 = SOM.transpose()*Hcore*SOM;
}

void HartreeFock::SetDensityMatrix(){
    // Builds the density matrix from the occupied MOs
    // By summing over all the occupied spatial MOs

    // Diagonlize Fock Matrix
    UTILS::Diagonlize(&F0, &e0, &C0);
    // Transform eigenvectors onto original non orthogonal AO basis
    C0 = SOM*C0;

    for (int i=0; i< NUM_ORB; i++) {
        for (int j=0; j< NUM_ORB; j++) {
            C0(i,0) = - C0(i,0);
            C0(i,4) = - C0(i,4);
            C0(i,6) = - C0(i,6);
        }
    }

    double M;
    for (int i = 0; i < NUM_ORB; i++){
        for (int j = 0; j < NUM_ORB; j++){
            M = 0;
            for(int m=0; m < NUM_OCC; m++) {
                M += C0(i,m) * C0(j,m);
            }
            D0(i,j) = M;
        }
    }
}

void HartreeFock::SaveDensity(){
  // Copies contents of D0 into prev_D0
  UTILS::CopyMatrix(&D0, &prev_D0);
}

void HartreeFock::SaveEnergy(){
    prev_etot = etot;
}

bool HartreeFock::EConverg(){
    // checks for convergence of the engery value
    delE = (prev_etot - etot);
    return (delE < tol_e);
}

bool HartreeFock::DensConverg(){
    // Checks for convergence of the density matrix
    rmsD = UTILS::RMS(&prev_D0, &D0);
    return rmsD < tol_dens;
}

void HartreeFock::SetEnergy() {
    eelec = 0;
    for (int i = 0; i < NUM_ORB; i++){
        for (int j = 0; j < NUM_ORB; j++){
            eelec += D0(i,j) * (Hcore(i,j) + F0(i,j));
        }
    }
    etot = eelec + enuc;
}

void HartreeFock::Iterate( int max_it ){

    int it = 0;
    cout << " -------------------------------------------- " << endl;
    cout << "Iter\t\t" << "E(elec)\t\t" << "E(tot)\t\t" << "Delta(E)\t\t" << "RMS(D)\t\t" << endl;
    cout << " -------------------------------------------- " << endl;

    while (not (EConverg() && DensConverg()) && it < max_it) {

        SaveDensity();
        SaveEnergy();

        SetFock();
        SetDensityMatrix();
        SetEnergy();

        if (it == 0) {

          cout << " -------------------------------------------- " << endl;
          cout << it << " Iteration Fock\t" << F0 << endl;
          cout << " -------------------------------------------- " << endl;
          cout << it << " Iteration Density\t" << D0 << endl;
          cout << " -------------------------------------------- " << endl;
        }

        cout << it << "\t\t" << eelec << "\t\t" << etot << "\t\t" << delE << "\t\t" << rmsD << endl;
        cout << " -------------------------------------------- " << endl;

        it ++;
    }
}

void HartreeFock::SetFock(){

    int ijkl, iklj;
    int ij, kl, ik, lj;

    for(int i =0; i < NUM_ORB; i++) {
        for(int j=0; j < NUM_ORB; j++) {

            F0(i,j) = Hcore(i,j);

            for(int k=0; k < NUM_ORB; k++) {
                for(int l=0; l < NUM_ORB; l++) {

                    ij = INDEX(i,j);
                    kl = INDEX(k,l);
                    ijkl = INDEX(ij,kl);
                    ik = INDEX(i,k);
                    lj = INDEX(l,j);
                    iklj = INDEX(ik,lj);

                    F0(i,j) += prev_D0(k,l) * ( TEI(ijkl) - 0.5 * TEI(iklj));
                }
            }
        }
    }
}
