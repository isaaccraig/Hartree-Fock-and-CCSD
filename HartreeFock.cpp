
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
using namespace Eigen;

HartreeFock::HartreeFock(string basisSet, double tol_dens, double tol_e){

    path = "data/" + basisSet + "/";
    READIN::val((path + "nBasis.dat").c_str(), &numBasisFunc);
    READIN::val((path + "nElectrons.dat").c_str(), &numElectrons);
    int temp = (numBasisFunc * (numBasisFunc + 1) / 2);
    int numMulliken = ((temp + 1) * temp)/2;

    S = new Eigen::MatrixXd(numBasisFunc, numBasisFunc);
    F0 = new Eigen::MatrixXd(numBasisFunc, numBasisFunc);
    V = new Eigen::MatrixXd(numBasisFunc, numBasisFunc);
    T = new Eigen::MatrixXd(numBasisFunc, numBasisFunc);
    Hcore = new Eigen::MatrixXd(numBasisFunc, numBasisFunc);

    errorVec = new Eigen::MatrixXd(numBasisFunc, numBasisFunc);
    SOM = new Eigen::MatrixXd(numBasisFunc, numBasisFunc);
    C0 = new Eigen::MatrixXd(numBasisFunc, numBasisFunc);
    FMO = new Eigen::MatrixXd(numBasisFunc, numBasisFunc);
    e0 = new Eigen::MatrixXd(numBasisFunc, numBasisFunc);
    D0 = new Eigen::MatrixXd(numBasisFunc, numBasisFunc);
    prev_D0 = new Eigen::MatrixXd(numBasisFunc, numBasisFunc);
    E = new Eigen::MatrixXd(numBasisFunc, 1);

    TEI_AO = new Eigen::MatrixXd(numMulliken, 1);
    TEI_MO = new Eigen::MatrixXd(numMulliken, 1);

    READIN::val((path + "enuc.dat").c_str(), &enuc);
    READIN::SymMatrix((path + "s.dat").c_str(),  S);
    READIN::SymMatrix((path + "t.dat").c_str(), T);
    READIN::SymMatrix((path + "v.dat").c_str(), V);
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

    //delete S;
    //delete V;
    //delete T;
    //delete Hcore;
    //delete E;
    //delete SOM;
    //delete F0;
    //delete FMO;
    //delete C0;
    //delete e0;
    //delete D0;
    //delete prev_D0;
    //delete TEI_MO;
    //delete TEI_AO;

    for (int i = 0; i < NUM_ERROR_MATRICES; i ++ ){
        //delete fockMatrices[i];
        //delete errorVectors[i];
    }

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
    cout << "Electronic Energy: \n" << eelec << endl;
}

bool HartreeFock::EConverg(){
    // checks for convergence of the engery value
    delE = (prev_etot - etot);
    return (delE < tol_e);
}

bool HartreeFock::DensConverg(){
    // Checks for convergence of the density matrix
    double val = 0;
    for (int i=0; i<numBasisFunc; i++){
        for (int j=0; j<numBasisFunc; j++){
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
    for (int i = 0; i < numBasisFunc; i++){
        for (int j = 0; j < numBasisFunc; j++){
            eelec += (*D0)(i,j) * ((*Hcore)(i,j) + (*F0)(i,j));
        }
    }
    etot = eelec + enuc;
}

void HartreeFock::SaveDensity(){
  for (int i = 0; i < numBasisFunc; i++) {
      for (int j = 0; j < numBasisFunc; j++) {
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

        cout << it << "\t\t" << setprecision(6) << eelec << "\t\t" << eelec + enuc << endl;

        it++;
    }
}

void HartreeFock::DIISIterate(){

    int it = 0;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Iter\t\t" << "Energy\t\t" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    while ( !(EConverg() && DensConverg()) ) {
        // Copy to check for convergence

        SaveEnergy();
        SaveDensity();
        Set_Fock();

        MatrixXd *F0copy = new MatrixXd(numBasisFunc, numBasisFunc);

        copyMatrix(F0, F0copy);
        fockMatrices[it%6] = F0copy;

        Set_Error(); // use normal fock in error calc
        errorVectors[it%6] = errorVec;
        if (it >= 2) {
            Extrapolate_Fock((it<NUM_ERROR_MATRICES)? it : NUM_ERROR_MATRICES);
        }
        Set_DensityMatrix(); // use extrapolated fock in density and energy calc
        Set_Energy();

        cout << it << "\t\t" << setprecision(6) << etot << endl;
        it ++;

    }
}

void HartreeFock::Extrapolate_Fock(int N) {

    MatrixXd *A = new MatrixXd(N + 1, N + 1);
    VectorXd *b = new VectorXd(N + 1);

    for (int i = 0; i < N + 1; i++){
        for (int j = 0; j <= i; j++) {
            if (i == N && j == N) {
                (*A)(i, j) = 0;
            } else if (i == N) {
                (*A)(i, j) = -1;
            } else {
                (*A)(i, j) = ((*errorVectors[i]) * (*errorVectors[j]).adjoint()).trace();
            }
            (*A)(j, i) = (*A)(i, j);
        }
        (*b)(i) = 0;
    }

    (*b)(N) = -1;
    VectorXd x = (*A).colPivHouseholderQr().solve(*b);
    setzero(F0);

    // last column solves for lagrange multiplier such that normalized, don't include
    // this is why bounds run to N, not N + 1 as they did prior
    for (int i = 0; i < N; i++) {
        (*F0) = x(i) * (*fockMatrices[i]);
    }

    //delete A;
    //delete b;
};

void HartreeFock::Set_Error(){
    // S :  AO-basis overlap matrix
    // D :  AO density used to make F (Fock)
    (*errorVec) = (*F0) * (*D0) * (*S) - (*S) * (*D0) * (*F0);
    (*errorVec) = (*SOM).adjoint() * (*errorVec) * (*SOM);

}

void HartreeFock::Set_Fock(){

    int ijkl, ikjl;
    int ij, kl, ik, jl;

    for(int i=0; i < numBasisFunc; i++) {
        for(int j=0; j < numBasisFunc; j++) {
            (*F0)(i,j) = (*Hcore)(i,j);
        }
    }

    for(int i=0; i < numBasisFunc; i++) {
        for(int j=0; j < numBasisFunc; j++) {

            for(int k=0; k < numBasisFunc; k++) {
                for(int l=0; l < numBasisFunc; l++) {

                    ij = compoundIndex(i,j);
                    kl = compoundIndex(k,l);
                    ijkl = compoundIndex(ij,kl);
                    ik = compoundIndex(i,k);
                    jl = compoundIndex(j,l);
                    ikjl = compoundIndex(ik,jl);

                    (*F0)(i,j) += (*prev_D0)(k,l) * (2.0 * (*TEI_AO)(ijkl) - (*TEI_AO)(ikjl));
                }
            }
        }
    }

    //cout << "------------------------------ Fock Matrix ------------------------------" << endl;
    //cout << (*F0) << endl;
    //(*F0) = (*SOM).transpose() * (*F0) * (*SOM);

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
    Diagonlize(F0, e0, C0);
    // Transform eigenvectors onto original AO basis
    (*C0) = (*SOM) * (*C0);

    setzero(D0);
    for (int i = 0; i < numBasisFunc; i++){
        for (int j = 0; j < numBasisFunc; j++){
            for(int m = 0; m < ((numElectrons + 4) * 0.5); m++) {
                (*D0)(i,j) += (*C0)(i,m) * (*C0)(j,m);
            }
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

void HartreeFock::DipoleMoment(){

    MatrixXd *mu_x= new MatrixXd(numBasisFunc, numBasisFunc);
    MatrixXd *mu_y= new MatrixXd(numBasisFunc, numBasisFunc);
    MatrixXd *mu_z= new MatrixXd(numBasisFunc, numBasisFunc);

    READIN::SymMatrix((path + "mux.dat").c_str(),mu_x);
    READIN::SymMatrix((path + "muy.dat").c_str(),mu_y);
    READIN::SymMatrix((path + "muz.dat").c_str(),mu_z);

    double mu_x_val;
    double mu_y_val;
    double mu_z_val;

    for (int mu = 0; mu < numBasisFunc; mu++){
        for (int nu = 0; nu < numBasisFunc; nu++){
            mu_x_val += (*D0)(mu,nu) * (*mu_x)(mu,nu);
            mu_y_val += (*D0)(mu,nu) * (*mu_y)(mu,nu);
            mu_z_val += (*D0)(mu,nu) * (*mu_z)(mu,nu);
        }
    }

    cout << "Mu-X = \t" << 2 * mu_x_val<< endl;
    cout << "Mu-Y = \t" << 2 * mu_y_val<< endl;
    cout << "Mu-Z = \t" << 2 * mu_z_val<< endl;
    cout << "Total Dipole Moment = \t" << pow(pow(2 * mu_x_val,2) + pow(2 * mu_y_val,2) + pow(2 * mu_z_val,2), 0.5) << endl;
    cout << "--------------------------------------------------------------------------------" << endl;


}

void HartreeFock::MullikenAnalysis(){

    // must know where each basis function is centered
    // this is not given !

    //Molecule mol("data/geom.txt", 0);
    int q = 0 ;
    int startOrb = 0;
    int endOrb = startOrb + 2;

    for (int i = 0; i < numBasisFunc; i ++) {
        //q = mol.zvals[1];
        for (int mu = startOrb; mu < endOrb; mu++) {
            q = q - 2 * ((*D0) * (*S))(mu, mu);
        }
        cout << "Charge on atom " << i << " = \n\n" << q << endl;
        startOrb = endOrb;
    }

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
    for (int i = 0; i< numBasisFunc; i++) {
        (*E)(i) = (*FMO)(i,i);
    }
}
