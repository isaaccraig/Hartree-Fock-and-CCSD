

#include "Read.hpp"

using namespace std;

void READIN::val(const char *filename, double *val) {
  FILE *input;
  input = fopen(filename, "r");
  fscanf(input, "%lf", val);
  fclose(input);
}

void READIN::SymMatrix(const char *filename, Matrix *M) {
  FILE *input;
  input = fopen(filename, "r");
  int i,j;
  double val;
  while (!feof(input)) {
    fscanf(input, "%d %d %lf", &i, &j, &val);
    (*M)(j-1,i-1) = val;
    (*M)(i-1,j-1) = val;
  }
  fclose(input);
}

void READIN::Mulliken(const char *filename, MullikenMatrix *TEI){
  FILE *input;
  input = fopen(filename, "r");
  int i, j, k, l, ij, kl, ijkl;
  double val;

  while(fscanf(input, "%d %d %d %d %lf", &i, &j, &k, &l, &val) != EOF) {

    ij = INDEX(i-1,j-1);
    kl = INDEX(k-1,l-1);
      
    ijkl = INDEX(ij,kl);

    if (ijkl > (*TEI).size()) {
        cout << "Error : Matrix Size Exceded" << endl;
        cout << (*TEI).size() << " = size" << endl;
        cout << ijkl << " = index" << endl;
        exit(-1);
    }
      
    (*TEI)(ijkl) = val;
      
  }
  fclose(input);
}

