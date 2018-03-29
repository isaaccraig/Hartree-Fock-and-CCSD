
#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include "Eigen/Dense"
#include "Eigen/Core"
#include "QuantumUtils.hpp"

struct Geometry {
    int natom;
    int *zvals;
    double **geom;
};

class READIN {

    public:
      static void val(const char *filename, double *val);
      static void val(const char *filename, int *val);
      static void SymMatrix(const char *filename, Eigen::MatrixXd *M);
      static void TEI(const char *filename, Eigen::MatrixXd *TEI);
      static void geometry(const char *filename, struct Geometry *g);

};
