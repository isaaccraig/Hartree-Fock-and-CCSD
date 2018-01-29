
#include "Utils.hpp"

using namespace std;

void UTILS::Diagonlize(Matrix *M, Matrix *evals, Matrix *evecs) {

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

void UTILS::CopyMatrix(Matrix *M1, Matrix *M2){
  for (int i = 0; i < NUM_ORB; i++) {
      for (int j = 0; j < NUM_ORB; j++) {
          (*M2)(i,j) = (*M1)(i,j);
      }
  }
}

void UTILS::SymmetricOrth(Matrix *S, Matrix *SOM) {
    // Diagonlizes S, such that S^-1/2 can be easily
    // calculated as L * U^-1/2 * L where L are the
    // evecs and U is the diagonal eigenvalue matrix

    Matrix eval;
    Matrix evec;
    Diagonlize(S, &eval, &evec);

    for(int i=0; i < NUM_ORB; i++) {
        eval(i,i) = pow(eval(i,i),-0.5);
    }

    (*SOM) = evec * eval * evec.transpose();
}

double UTILS::RMS(Matrix *M0, Matrix *M1){
    // Returns RMS difference between two matricies
    double val = 0;
    for (int i=0; i < NUM_ORB; i++){
        for (int j=0; j < NUM_ORB; j++){
            val += pow((*M0)(i,j) - (*M1)(i,j), 2);
        }
    }
    return pow(val, 0.5);
}

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

  for (int i = 0; i < (*M).size(); i++) {
    for (int j = 0; j < (*M).size(); j++) {
      (*M)(i,j) = FILLER;
    }
  }

  while (!feof(input)) {
    fscanf(input, "%d %d %lf", &i, &j, &val);
    (*M)(j-1,i-1) = val;
    (*M)(i-1,j-1) = val;
  }
  fclose(input);
}

void READIN::TEI(const char *filename, TEIMatrix *TEI){
  FILE *input;
  input = fopen(filename, "r");
  int i, j, k, l, ij, kl, ji, lk;
  double val;

  for (int i = 0; i < (*TEI).size(); i++) {
      (*TEI)(i) = FILLER;
  }

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
