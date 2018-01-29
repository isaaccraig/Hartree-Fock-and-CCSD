
#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include "Eigen/Dense"
#include "Eigen/Core"

#define NUM_ORB 7
#define NUM_ELEC 10
#define NUM_OCC NUM_ELEC/2
#define NUM_MUL 407
//NUM_ORB*(NUM_ORB+1)/2

inline int INDEX(int i, int j) {
  return ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)));
}

typedef Eigen::Matrix<double, NUM_ORB, NUM_ORB> Matrix;
typedef Eigen::Matrix<double, NUM_MUL, 1> TEIMatrix;
typedef Eigen::Matrix<double, NUM_ORB, 1> Vector;
typedef Eigen::Matrix<double, NUM_MUL, NUM_MUL> Double_Matrix;

void setzero(Matrix *M);
void setzero(TEIMatrix *M);
void setzero(Double_Matrix *M);

class READIN {

    public:
      static void val(const char *filename, double *val);
      static void SymMatrix(const char *filename, Matrix *M);
      static void TEI(const char *filename, TEIMatrix *TEI);

};
