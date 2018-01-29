
//
//  main.cpp
//  HartreeFock
//
//  Created by Isabel Craig on 1/26/18.
//  Copyright © 2018 Isabel Craig. All rights reserved.
//

#include "HartreeFock.hpp"
#include "Read.hpp"

using namespace std;

HartreeFock::HartreeFock(double tol_dens, double tol_e){
    
    READIN::val("data/enuc.dat", &enuc);
    READIN::SymMatrix("data/overlap.dat",  &S);
    READIN::SymMatrix("data/kinetic.dat", &T);
    READIN::SymMatrix("data/anuc.dat", &V);
    READIN::Mulliken("data/eri.dat", &TEI);
    
    Hcore = T + V;
    
    SymmetricOrth();              // Symmetric Orthogalization Matrix
    Set_InitialFock();                // Build Initial Guess Fock Matrix
    Set_DensityMatrix();              // Build Initial Density Matrix using occupied MOs
    Set_Energy();                  // Compute the Initial SCF Energy
        
}

void Diagonlize(Matrix *M, Matrix *evals, Matrix *evecs) {
    
    Eigen::SelfAdjointEigenSolver<Matrix> solver(*M);
    *evecs = solver.eigenvectors();
    Vector evals_vec = solver.eigenvalues();
    for (int i=0; i<NUM_ORB; i++) {
        for (int j=0; j<NUM_ORB; j++) {
            if (i==j) {
                (*evals)(i,i) = evals_vec(i);
            } else {
                (*evals)(i,j) = 0;
            }
        }
    }
}

void HartreeFock::print_state() {
    
    cout << "Nuclear repulsion energy = \n" << enuc << endl;
    cout << "Overlap Integrals: \n" << S << endl;
    cout << "Kinetic-Energy Integrals: \n" << T << endl;
    cout << "Nuclear Attraction Integrals: \n" << V << endl;
    cout << "Core Hamiltonian: \n" << Hcore << endl;
    cout << "S^-1/2 Matrix: \n" << SOM << endl;
    cout << "Initial F' Matrix: \n" << F0 << endl;
    cout << "Initial C Matrix: \n" << C0 << endl;
    cout << "Initial Density Matrix: \n" << D0 << endl;
    cout << "Initial Energy: \n" << etot << endl;
}

bool HartreeFock::EConverg(){
    // checks for convergence of the engery value
    delE = ( (prev_etot > etot) ? (prev_etot - etot) : (etot - prev_etot) );
    return (delE < tol_e);
}

bool HartreeFock::DensConverg(){
    // Checks for convergence of the density matrix
    double val = 0;
    for (int i=0; i<NUM_ORB; i++){
        for (int j=0; j<NUM_ORB; j++){
            val += pow(prev_D0(i,j) - D0(i,j), 2);
        }
    }
    rmsD = pow(val, 0.5);
    return (rmsD < tol_dens);
}

void HartreeFock::Set_Energy() {
    // Sum over all atomic orbitals
    // of DensityMatrix * (Hcore + Fock)
    eelec = 0;
    for (int i = 0; i < NUM_ORB; i++){
        for (int j = 0; j < NUM_ORB; j++){
            eelec += D0(i,j) * (Hcore(i,j) + F0(i,j));
        }
    }
    etot = eelec + enuc;
}

void HartreeFock::Iterate( int maxit ){
    
    int it = 0;
    cout << "Iteration\t\t" << "E(elec)\t\t" << "E(tot)\t\t" << "Delta(E)\t\t" << "RMS(D)\t\t" << endl;
    while ( (not EConverg()) && (not DensConverg())) {
        // Copy to check for convergence
        for (int i = 0; i < NUM_ORB; i++) {
            for (int j = 0; j < NUM_ORB; j++) {
                prev_D0(i,j) = D0(i,j);
            }
        }
        
        prev_etot = etot;
        Set_Fock();
        Set_DensityMatrix();
        Set_Energy();
        
        cout << it << "\t\t" << eelec << "\t\t" << etot << "\t\t" << delE << "\t\t" << rmsD << endl;
        
        it ++;
        if (it < maxit) {
            exit(-1);
        }
    }

}

void HartreeFock::Set_Fock(){
    
    int ijkl, ikjl;
    int ij, kl, ik, jl;
    
    for(int i=0; i < NUM_ORB; i++) {
        for(int j=0; j < NUM_ORB; j++) {
            
            F0(i,j) = Hcore(i,j);
            
            for(int k=0; k < NUM_ORB; k++) {
                for(int l=0; l < NUM_ORB; l++) {
                    
                    ij = INDEX(i,j);
                    kl = INDEX(k,l);
                    ijkl = INDEX(ij,kl);
                    ik = INDEX(i,k);
                    jl = INDEX(j,l);
                    ikjl = INDEX(ik,jl);
                    
                    F0(i,j) += D0(k,l) * (2.0 * TEI[ijkl] - TEI[ikjl]);
                }
            }
        }
    }
}

void HartreeFock::SymmetricOrth() {
    // Diagonlizes S, such that S^-1/2 can be easily
    // calculated as L * U^-1/2 * L where L are the
    // evecs and U is the diagonal eigenvalue matrix
    
    Matrix eval;
    Matrix evec;
    Diagonlize(&S, &eval, &evec);
    
    for(int i=0; i < NUM_ORB; i++) {
        eval(i,i) = pow(eval(i,i),-0.5);
    }
    
    SOM = evec * eval * evec.transpose();
}

void HartreeFock::Set_InitialFock(){
    // forms an intial guess fock matrix in orthonormal AO using
    // the core hamiltonian as a Guess, such that
    // Fock = transpose (S^-1/2) * Hcore * S^-1/2
    F0 = SOM.transpose()*Hcore*SOM;
}

void HartreeFock::Set_DensityMatrix(){
    // Builds the density matrix from the occupied MOs
    // By summing over all the occupied spatial MOs
    
    Diagonlize(&F0, &e0, &C0);      // Diagonlize Fock Matrix
    C0 = SOM*C0;                    // Transform eigenvectors onto original non orthogonal AO basis
    
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

