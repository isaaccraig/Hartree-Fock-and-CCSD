

#include "Read.hpp"

using namespace std;

void setzero(Matrix *M){
  for (int i=0; i<M->rows(); i++){
      for (int j=0; j<M->cols(); j++){
        (*M)(i,j) = 0;
      }
  }
}

void setzero(TEIMatrix *M){
  for (int i=0; i<M->rows(); i++){
      for (int j=0; j<M->cols(); j++){
        (*M)(i,j) = 0;
      }
  }
}

void setzero(Double_Matrix *M){
  for (int i=0; i<M->rows(); i++){
      for (int j=0; j<M->cols(); j++){
        (*M)(i,j) = 0;
      }
  }
}

void READIN::val(const char *filename, double *val) {
  FILE *input;
  input = fopen(filename, "r");
  fscanf(input, "%lf", val);
  fclose(input);
}

void READIN::SymMatrix(const char *filename, Matrix *M) {
  FILE *input;
  setzero(M);
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

void READIN::TEI(const char *filename, TEIMatrix *TEI){
  FILE *input;
  setzero(TEI);
  input = fopen(filename, "r");
  int i, j, k, l, ij, kl, ji, lk;
  double val;

  while(fscanf(input, "%d %d %d %d %lf", &i, &j, &k, &l, &val) != EOF) {

    ij = INDEX(i-1,j-1);
    kl = INDEX(k-1,l-1);
    ji = INDEX(j-1,i-1);
    lk = INDEX(l-1,k-1);

    (*TEI)(INDEX(ij,kl)) = val;
    (*TEI)(INDEX(ji,kl)) = val;
    (*TEI)(INDEX(ij,lk)) = val;
    (*TEI)(INDEX(ji,lk)) = val;
    (*TEI)(INDEX(kl,ij)) = val;
    (*TEI)(INDEX(kl,ji)) = val;
    (*TEI)(INDEX(lk,ji)) = val;
    (*TEI)(INDEX(lk,ij)) = val;

  }
  fclose(input);
}
