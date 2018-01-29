
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
#define NUM_MUL 406

inline int INDEX(int i, int j) {
  return ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)));
}

typedef Eigen::Matrix<double, NUM_ORB, NUM_ORB> Matrix;
typedef Eigen::Matrix<double, NUM_MUL, 1> MullikenMatrix;
typedef Eigen::Matrix<double, NUM_ORB, 1> Vector;

class READIN {

    public:
      static void val(const char *filename, double *val);
      static void SymMatrix(const char *filename, Matrix *M);
      static void Mulliken(const char *filename, MullikenMatrix *TEI);
    
    private:
      static void Test_Mulliken(MullikenMatrix *TEI);

};
