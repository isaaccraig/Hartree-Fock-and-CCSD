
#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <set>
#include <exception>
#include "Eigen/Dense"
#include "Eigen/Core"

#define NUM_ORB 14
#define NUM_ELEC 10
#define NUM_OCC (NUM_ELEC/2)
#define NUM_MUL 477

class SanityCheckException: public std::exception {
  std::string msg;
  public:
    SanityCheckException(const char* msg);
    virtual const char* what() const throw() { return (this->msg).c_str(); };
};

inline int compoundIndex(int i, int j) {
    return ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)));
}


void printMatrix(Eigen::MatrixXd *M, std::string name);
double**** allocate4DMatrix(int n);
void delete4DMatrix(double**** M, int n);
double getExchangeIntegral(Eigen::MatrixXd *TEI, int i, int j, int a, int b);
double getExchangeIntegral(double ****TEI, int i, int j, int a, int b);
void spinOrbitalFock(Eigen::MatrixXd *FSO, double**** TEI_MOspin, Eigen::MatrixXd *h, int n);
void Diagonlize(Eigen::MatrixXd *M, Eigen::MatrixXd *evals, Eigen::MatrixXd *evecs);
void Diagonlize(Eigen::MatrixXd *M, Eigen::MatrixXd *evals, Eigen::MatrixXd *evecs, Eigen::MatrixXd *S);
bool isNormalized(Eigen::MatrixXd *M, Eigen::MatrixXd *S);
void NormalizeColumns(Eigen::MatrixXd *M, Eigen::MatrixXd *S);
void SymmetricOrthMatrix(Eigen::MatrixXd *SOM, Eigen::MatrixXd *S);
void setzero(Eigen::MatrixXd *M);
void atomicToMolecularN8(Eigen::MatrixXd *TEI_MO, Eigen::MatrixXd *TEI_AO, Eigen::MatrixXd *C0);
void atomicToMolecularN5(Eigen::MatrixXd *TEI_MO, Eigen::MatrixXd *TEI_AO, Eigen::MatrixXd *C0);
void molecularToMolecularSpin(double**** TEI_MOspin, Eigen::MatrixXd *TEI_MO, int n);
double MP2_Energy(Eigen::MatrixXd *TEI_MO, Eigen::MatrixXd *E, int nElectrons);
void copyMatrix(Eigen::MatrixXd *M, Eigen::MatrixXd *N);
