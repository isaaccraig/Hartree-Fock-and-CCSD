
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <string>

#include "molecule.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

namespace utils { // Error and Warning Calls
  /* base Error to give message and halt */
  inline void error( const char *msg ) {
    fprintf(stderr, "Error : %s\n", msg);
    exit(-1);
  }
  /* Warning to give message and continue */
  inline void warning( const char *msg ) {
    fprintf(stderr, "Warning : %s\n", msg);
  }

  inline Eigen::Vector2i other_indeces(int index){
    Eigen::Vector2i indeces;
    int end = 0;
    for (int i=0; i<3; i++){
      if (not (i == index)){
        indeces(end) = i;
        end++;
      }
    }
    return indeces;
  }

};

void Molecule::print_geom(){
  for(int i=0; i < natom; i++)
    printf("%d %8.5f %8.5f %8.5f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
}

void Molecule::translate(double x, double y, double z){
  for(int i=0; i < natom; i++) {
     geom[i][0] += x;
     geom[i][1] += y;
     geom[i][2] += z;
  }
}

double Molecule::bond(int a, int b) {
  return sqrt( (geom[a][0]-geom[b][0])*(geom[a][0]-geom[b][0])
             + (geom[a][1]-geom[b][1])*(geom[a][1]-geom[b][1])
             + (geom[a][2]-geom[b][2])*(geom[a][2]-geom[b][2]) );
}

// Returns the value of the unit vector between atoms a and b
// in the cart direction (cart=0=x, cart=1=y, cart=2=z)
double Molecule::unit(int cart, int a, int b){
  return - (geom[a][cart] - geom[b][cart])/bond(a,b);
}

void Molecule::unitvec(int a, int b, Eigen::Vector3d v){
  v(0) = unit(0, a, b);
  v(1) = unit(1, a, b);
  v(2) = unit(2, a, b);
}

// Returns the angle between atoms a, b, and c in radians
double Molecule::angle(int a, int b, int c) {
  return acos(unit(0,b,a) * unit(0,b,c) + unit(1,b,a) * unit(1,b,c) + unit(2,b,a) * unit(2,b,c));
}

double Molecule::outofplane(int a, int b, int c, int d){
  double theta = ((unit(1,c,b) * unit(2,c,d)) - (unit(2,c,b) * unit(1,c,d)))/sin(angle(b,c,d)) * unit(0,c,a) -
                ( (unit(0,c,b) * unit(2,c,d)) - (unit(2,c,b) * unit(0,c,d)))/sin(angle(b,c,d)) * unit(1,c,a) +
                ( (unit(0,c,b) * unit(1,c,d)) - (unit(1,c,b) * unit(0,c,d)))/sin(angle(b,c,d)) * unit(2,c,a);
  if(theta < -1.0) return asin(-1.0);
  else if(theta > 1.0) return asin(1.0);
  else return asin(theta);
}

double Molecule::torsional(int i, int j, int k, int l){

  Eigen::Vector3d u_ij;
  unitvec(i,j,u_ij);

  Eigen::Vector3d u_jk;
  unitvec(j,k, u_jk);

  Eigen::Vector3d u_kl;
  unitvec(k,l, u_kl);

  Eigen::Vector3d t1;
  t1 = u_ij.cross(u_jk) / sin(angle(i,j,k));

  Eigen::Vector3d t2;
  t2 = u_jk.cross(u_kl) / sin(angle(j,k,l));

  return acos(t1.dot(t2));
}

double Molecule::masscenter(int cart){

  double cm = 0;
  double M = 0;

  for(int i=0; i < natom; i++) M += atomic_mass[(int) zvals[i]];

  for (int i=0; i<natom; i++){
    cm += atomic_mass[(int) zvals[i]] * geom[i][cart];
  }
  return cm/M;
}

void Molecule::inertia(Eigen::Matrix3d *evals, Eigen::Matrix3d *evecs){

  Eigen::Matrix3d I;
  Eigen::Vector2i ind;

  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      if (i==j){
        for (int n=0; n<natom; n++) {
          ind = utils::other_indeces(i);
          I(i,j) = atomic_mass[(int) zvals[i]] * (pow(geom[n][ind(0)], 2) + pow(geom[n][ind(1)], 2));
        }
      }
      else if(j>i) {
        I(j,i) = I(i,j);
      }
      else{
        for (int n=0; n<natom; n++) {
          I(i,j) = atomic_mass[(int) zvals[i]] * geom[n][i] * geom[n][j];
        }
      }
    }
  }

  cout << "Moment of inertia tensor (amu bohr^2):" << endl;
  cout << I << endl;
  // Diagonlize using Eigen

  Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
  *evecs = solver.eigenvectors();
  *evals = solver.eigenvalues();

  cout << "Principle Moments of Interia:" << endl;
  cout << evals << endl;

  cout << "\nPrincipal moments of inertia (amu * AA^2):\n";
  cout << evals * 0.529177249 * 0.529177249 << endl;

  cout << "\nPrincipal moments of inertia (g * cm^2):\n";
  cout << evals * 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8 << endl;

  cout << "Rotational Constants:" << endl;
  double coef = h / (8 * 3.14159 * 3.14159 * c);
  cout << coef * evals(0) << coef * evals(1) << coef * evals(2) << endl;


}

void Molecule::set_geom(const char *filename){
  // open filename
  ifstream input(filename);
  if (not input.is_open()) {
    std::string message = "Unable to open file : ";
    message += filename;
    utils::warning(message.c_str());
  }
  // read the number of atoms from filename
  input >> natom;
  zvals = new int[natom];

  // Dynamically Allocate Space for geom
  geom = new double* [natom];
  for(int i=0; i < natom; i++)
    geom[i] = new double[3];

  // read in the geometry from filename
  for(int i=0; i < natom; i++)
    input >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];

  input.close();
}

Molecule::Molecule(const char *filename, int q){
  charge = q;
  atomic_mass = {
    0.0000000,
    1.00782503223,
    4.00260325413,
    6.0151228874,
    9.012183065,
    11.00930536,
    12.00000000000,
    14.00307400443,
    15.99491461957,
    18.99840316273,
    19.9924401762,
    22.9897692820,
    23.985041697,
    26.98153853,
    27.97692653465,
    30.97376199842,
    31.9720711744,
    34.968852682,
    39.9623831237,
    38.9637064864
    };
  // open filename
    set_geom(filename);
}

Molecule::~Molecule(){
  delete[] zvals;

  for(int i=0; i < natom; i++)
    delete[] geom[i];
  delete[] geom;
}

