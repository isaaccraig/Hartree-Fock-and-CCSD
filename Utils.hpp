
#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cassert>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

#define NUM_ORB 7
#define NUM_TRI NUM_ORB*(NUM_ORB+1)/2
#define NUM_ELEC 10
#define NUM_OCC NUM_ELEC/2
#define NUM_MUL NUM_TRI*(NUM_TRI+1)/2
#define FILLER 0

inline int INDEX(int i, int j) {
  return ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)));
}

typedef Eigen::Matrix<double, NUM_ORB, NUM_ORB> Matrix;
typedef Eigen::Matrix<double, NUM_MUL, 1> TEIMatrix;
typedef Eigen::Matrix<double, NUM_ORB, 1> Vector;

class UTILS {

    public:
      static void Diagonlize(Matrix *M, Matrix *evals, Matrix *evecs);
      static void SymmetricOrth(Matrix *S, Matrix *SOM);
      static double RMS(Matrix *M0, Matrix *M1);
      static void CopyMatrix(Matrix *M0, Matrix *M1);
};

class READIN {

    public:
      static void val(const char *filename, double *val);
      static void SymMatrix(const char *filename, Matrix *M);
      static void TEI(const char *filename, TEIMatrix *TEI);

};
