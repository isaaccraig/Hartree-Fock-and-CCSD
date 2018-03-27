
//
//  main.cpp
//  HartreeFock
//
//  Created by Isabel Craig on 1/26/18.
//  Copyright © 2018 Isabel Craig. All rights reserved.
//

#include "HartreeFock.hpp"
#include "QuantumUtils.hpp"
#include "Read.hpp"
#include <cmath>

using namespace std;

HartreeFock::HartreeFock(double tol_dens, double tol_e){

    string path = "data/STO3G_Water/";

    S = new Eigen::MatrixXd(NUM_ORB, NUM_ORB);
    F0 = new Eigen::MatrixXd(NUM_ORB, NUM_ORB);
    V = new Eigen::MatrixXd(NUM_ORB, NUM_ORB);
    T = new Eigen::MatrixXd(NUM_ORB, NUM_ORB);

    Hcore = new Eigen::MatrixXd(NUM_ORB, NUM_ORB);
    SOM = new Eigen::MatrixXd(NUM_ORB, NUM_ORB);
    C0 = new Eigen::MatrixXd(NUM_ORB, NUM_ORB);
    FMO = new Eigen::MatrixXd(NUM_ORB, NUM_ORB);

    e0 = new Eigen::MatrixXd(NUM_ORB, NUM_ORB);
    D0 = new Eigen::MatrixXd(NUM_ORB, NUM_ORB);
    prev_D0 = new Eigen::MatrixXd(NUM_ORB, NUM_ORB);

    E = new Eigen::MatrixXd(NUM_ORB, 1);

    TEI_AO = new Eigen::MatrixXd(NUM_MUL, 1);
    TEI_MO = new Eigen::MatrixXd(NUM_MUL, 1);

    READIN::val((path + "enuc.dat").c_str(), &enuc);
    READIN::SymMatrix((path + "overlap.dat").c_str(),  S);
    READIN::SymMatrix((path + "kinetic.dat").c_str(), T);
    READIN::SymMatrix((path + "anuc.dat").c_str(), V);
    READIN::TEI((path + "eri.dat").c_str(), TEI_AO);

    this->tol_dens = tol_dens;
    this->tol_e = tol_e;

    *Hcore = *T + *V;

    SymmetricOrthMatrix(SOM, S);                // Symmetric Orthogalization Matrix
    Set_InitialFock();                          // Build Initial Guess Fock Matrix
    Set_DensityMatrix();                        // Build Initial Density Matrix using occupied MOs
    Set_Energy();                               // Compute the Initial SCF Energy

}

HartreeFock::~HartreeFock() {

    delete S;
    delete V;
    delete T;
    delete Hcore;
    delete E;
    delete SOM;
    delete F0;
    delete FMO;
    delete C0;
    delete e0;
    delete D0;
    delete prev_D0;

    delete TEI_MO;
    delete TEI_AO;
}

void HartreeFock::print_state() {

    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "------------------------ Hartree Fock w/ MP2 Correction ------------------------" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Nuclear repulsion energy = " << enuc << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Overlap Integrals: \n" << (*S) << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Kinetic-Energy Integrals: \n" << (*T) << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Nuclear Attraction Integrals: \n" << (*V) << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Core Hamiltonian: \n" << (*Hcore) << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Symmetric Orthogalization Matrix: \n" << (*SOM) << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Fock Matrix (In Orthogalized Basis): \n" << (*F0) << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "MO Coefficent Matrix: \n" << (*C0) << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Density Matrix: \n" << (*D0) << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Energy: \n" << etot << endl;
}

bool HartreeFock::EConverg(){
    // checks for convergence of the engery value
    delE = (prev_etot - etot);
    return (delE < tol_e);
}

bool HartreeFock::DensConverg(){
    // Checks for convergence of the density matrix
    double val = 0;
    for (int i=0; i<NUM_ORB; i++){
        for (int j=0; j<NUM_ORB; j++){
            val += pow((*prev_D0)(i,j) - (*D0)(i,j), 2);
        }
    }
    rmsD = pow(val, 0.5);
    return (rmsD < tol_dens);
}

void HartreeFock::CheckEnergy(){

  double expected = -74.991229564312;
  double percent_off = 100 * (EMP2 + etot - expected)/expected;
  cout << "--------------------------------------------------------------------------------" << endl;
  cout << percent_off << " percent off from expected results" << endl;
  cout << (EMP2 + etot - expected) << " off from expected results" << endl;
  cout << "--------------------------------------------------------------------------------" << endl;

}

void HartreeFock::Set_Energy() {
    // Sum over all atomic orbitals
    // of DensityMatrix * (Hcore + Fock)
    eelec = 0;
    for (int i = 0; i < NUM_ORB; i++){
        for (int j = 0; j < NUM_ORB; j++){
            eelec += (*D0)(i,j) * ((*Hcore)(i,j) + (*F0)(i,j));
        }
    }
    etot = eelec + enuc;
}

void HartreeFock::SaveDensity(){
  for (int i = 0; i < NUM_ORB; i++) {
      for (int j = 0; j < NUM_ORB; j++) {
          (*prev_D0)(i,j) = (*D0)(i,j);
      }
  }
}

void HartreeFock::SaveEnergy(){
    prev_etot = etot;
}

void HartreeFock::Iterate(){

    int it = 0;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Iter\t\t" << "Energy\t\t" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    while ( !(EConverg() && DensConverg()) ) {
        // Copy to check for convergence

        SaveEnergy();
        SaveDensity();
        Set_Fock();
        Set_DensityMatrix();
        Set_Energy();

        cout << it << "\t\t" << setprecision(6) << etot << endl;

        it ++;
    }
}

void HartreeFock::Set_Fock(){

    int ijkl, ikjl;
    int ij, kl, ik, jl;

    for(int i=0; i < NUM_ORB; i++) {
        for(int j=0; j < NUM_ORB; j++) {

            (*F0)(i,j) = (*Hcore)(i,j);

            for(int k=0; k < NUM_ORB; k++) {
                for(int l=0; l < NUM_ORB; l++) {

                    ij = compoundIndex(i,j);
                    kl = compoundIndex(k,l);
                    ijkl = compoundIndex(ij,kl);
                    ik = compoundIndex(i,k);
                    jl = compoundIndex(j,l);
                    ikjl = compoundIndex(ik,jl);

                    (*F0)(i,j) += (*D0)(k,l) * (2.0 * (*TEI_AO)(ijkl) - (*TEI_AO)(ikjl));
                }
            }
        }
    }
}

void HartreeFock::Set_InitialFock(){
    // forms an intial guess fock matrix in orthonormal AO using
    // the core hamiltonian as a Guess, such that
    // Fock = transpose (S^-1/2) * Hcore * S^-1/2
    (*F0) = (*SOM).transpose() * (*Hcore) * (*SOM);
}

void HartreeFock::Set_DensityMatrix(){
    // Builds the density matrix from the occupied MOs
    // By summing over all the occupied spatial MOs

    // Diagonlize Fock Matrix
    Diagonlize(F0, e0, C0);

    // Transform eigenvectors onto original AO basis
    (*C0) = (*SOM) * (*C0);

    double M;
    for (int i = 0; i < NUM_ORB; i++){
        for (int j = 0; j < NUM_ORB; j++){
            M = 0;
            for(int m=0; m < NUM_OCC; m++) {
                M += (*C0)(i,m) * (*C0)(j,m);
            }
            (*D0)(i,j) = M;
        }
    }
}

void HartreeFock::Set_MOBasisFock() {

    // Tests that the resultant Fock matrix is diagonal in the MO basis
    // orbital elements should be diagonal elements since Fi |xi> = ei |xi>
    // therefore Fij = <xi|F|xj> = ei * dij

    // Convert from AO to MO using LCAO-MO coefficents
    //      MO(i)   = Sum over v of C(m,i) * AO(m)
    // Fij = Sum over m,v of C(m,j) * C(v,i) <psi m|F|psi v>
    //              = Sum over m,v of C(m,j) * C(v,i) * F(m,v)

    //cout << SOM.transpose()*F0*SOM << endl;

    setzero(FMO);
    (*FMO) = (*C0).transpose() * (*F0) * (*C0);
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "MO Basis Fock Matrix:\n" << (*FMO) << endl;

}

void HartreeFock::MP2_Correction(){

    Set_MOBasisFock();
    Set_OrbitalEnergy();
    atomicToMolecularN8(TEI_MO, TEI_AO, C0);
    EMP2 = MP2_Energy(TEI_MO, E);
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "MP2 Correction Energy :" << EMP2 << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Corrected Energy :" << EMP2 + etot << endl;
    cout << "--------------------------------------------------------------------------------" << endl;

}

void HartreeFock::Set_OrbitalEnergy(){
    // Diagonal Elements of The Fock Operator
    // in the MO Bais are the orbital Energy values
    for (int i = 0; i< NUM_ORB; i++) {
        (*E)(i) = (*FMO)(i,i);
    }
}

