
#include "QuantumUtils.hpp"
#include <stdexcept>
#include <stdio.h>
using namespace std;
using namespace Eigen;


double**** allocate4DMatrix(int n) {
    double**** W = new double*** [n];
    for(int i=0; i < n; i++) {
        W[i] = new double** [n];
        for(int j=0; j < n; j++) {
            W[i][j] = new double* [n];
            for(int k=0; k < n; k++) {
                W[i][j][k] = new double[n];
            }
        }
    }
    for(int i=0; i < n; i++) {
        for(int j=0; j < n; j++) {
            for(int k=0; k < n; k++) {
                for(int l=0; l < n; l++) {
                    W[i][j][k][l] = 0;
                }
            }
        }
    }
    return W;
}

void delete4DMatrix(double**** W, int n) {
    for(int i=0; i < n; i++) {
        for(int j=0; j < n; j++) {
            for(int k=0; k < n; k++) {
                delete W[i][j][k];
            }
            delete W[i][j];
        }
        delete W[i];
    }
    delete W;
}

/** Diagonlizes S, such that S^-1/2 can be easily
calculated as L * U^-1/2 * L where L are the
evecs and U is the diagonal eigenvalue matrix
**/
void SymmetricOrthMatrix(MatrixXd *SOM, MatrixXd *S) {

    if (S == NULL || SOM == NULL ){
        throw invalid_argument( "received null pointer in SymmetricOrthMatrix" );
    }

    int numOrb = (*S).rows();
    MatrixXd *eval = new MatrixXd(numOrb, numOrb);
    MatrixXd *evec = new MatrixXd(numOrb, numOrb);

    Diagonlize(S, eval, evec);
    for(int i=0; i < numOrb; i++) {
        (*eval)(i,i) = pow((*eval)(i,i),-0.5);
    }
    *SOM = (*evec) * (*eval) * (*evec).transpose();

    delete eval;
    delete evec;
}


/**
Diagonalizes M, returning a matrix of right eigenvectors
and a diagonal of eigenvalues
**/
void Diagonlize(MatrixXd *M, MatrixXd *evals, MatrixXd *evecs) {
    if (evals == NULL || evecs == NULL || M == NULL ){
        throw invalid_argument( "received null pointer in Diagonlize" );
    }
    int numOrb = (*M).rows();
    SelfAdjointEigenSolver<MatrixXd> solver(*M);
    *evecs = solver.eigenvectors();
    auto evals_vec = solver.eigenvalues();
    for (int i=0; i<numOrb; i++) {
        for (int j=0; j<numOrb; j++) {
            if (i==j) {
                (*evals)(i,i) = evals_vec(i);
            } else {
                (*evals)(i,j) = 0;
            }
        }
    }
}

/**
Diagonalizes M, returning a matrix of right eigenvectors
and a diagonal of eigenvalues
**/
void Diagonlize(MatrixXd *M, MatrixXd *evals, MatrixXd *evecs, MatrixXd *S) {
    int numOrb = (*M).rows();
    SelfAdjointEigenSolver<MatrixXd> solver(*M);
    *evecs = solver.eigenvectors();
    auto evals_vec = solver.eigenvalues();
    for (int i=0; i<numOrb; i++) {
        for (int j=0; j<numOrb; j++) {
            if (i==j) {
                (*evals)(i,i) = evals_vec(i);
            } else {
                (*evals)(i,j) = 0;
            }
        }
    }
    NormalizeColumns(evecs, S);
}

/**
Normalizes the Columns, c, of Matrix M such that
c' * c = 1
**/
void NormalizeColumns(MatrixXd *M, MatrixXd *S){
  // Assert that c * S * c = 1
    int numOrb = (*M).rows();
    VectorXd *eigenvec = new VectorXd(numOrb);
    double norm_val;
    for (int i = 0; i < (*M).cols(); i++){
        // pull out eivenvector
        (*eigenvec) = (*M).col(i);
        // normilization value is c * S * c
        norm_val = (*eigenvec).transpose() * (*S) * (*eigenvec);
        // eigenvector is now 1/root(normilization value) * eigenvector
        (*eigenvec) = (*eigenvec) * pow(norm_val, -0.5);
        (*M).col(i) = (*eigenvec);

    }
    if (!isNormalized(M,S)){
        throw new SanityCheckException("Unsuccessful Normalization");
    }
}

/**
returns true if for all columns c in M,
c' * S * c = 1, where S is the overlap matrix
**/
bool isNormalized(MatrixXd *M, MatrixXd *S){
    auto prod = (*M).transpose() * (*S) * (*M);
    auto diag = (prod).diagonal();
    for (int i = 0; i < (diag).rows(); i++){
          if ((diag)(i) != 1) {
             return false;
          }
    }
    return true;
}

void setzero(MatrixXd *M){
    if (M == NULL){
        throw invalid_argument( "received null pointer in setZero" );
    }
  for (int i=0; i<M->rows(); i++){
      for (int j=0; j<M->cols(); j++){
        (*M)(i,j) = 0;
      }
  }
}

/**
custom exception to be thrown when intermediate Sanity
checks fail, constructed with message;
**/

SanityCheckException::SanityCheckException(const char* msg) {
    string errMsg = "Sanity Check Failure :";
    this->msg = errMsg + msg;
}


/**
Translates the TEI from the MO basis to the
Molecular Spin Orbital Basis

@NOTE:
  requires all occupied orbitals before virtual orbitals
  alternating between spins

runtime : O(N^4)

            (p%2 == r%2) and (q%2 == s%2) factors checks
            if the spin of the two AO wavefunctions for
            each cordinate are equal, as if they are not,
            the < pq | rs > value would be zero.

            @notation : <pq||rs> = <pq|rs> - <ps|qr> = spin orbital TEI (MO basis)
            @notation : (pq|rs) = orbital TEI (MO basis)

            <pq|rs> = (pq|rs) <spin_p|spin_r> <spin_q|spin_s>
                    = TEI_MO(pqrs) * (p%2 == r%2) * (q%2 == s%2)

            <pq||rs> =  <pq|rs> - <ps|qr>
                     =  TEI_MO(pqrs) * (p%2 == r%2) * (q%2 == s%2) -
                        TEI_MO(psqr) * (p%2 == s%2) * (r%2 == q%2)

**/
void molecularToMolecularSpin(double ****TEI_MOspin, MatrixXd *TEI_MO, int numOrb) {

    if (TEI_MO == NULL || TEI_MOspin == NULL){
        throw invalid_argument( "received null pointer in molecularToMolecularSpin" );
    }
    /**
        assumes
        even-numbered orbitals are spin up
        odd-numbered are spin down
    **/

    for (int p=0; p < numOrb; p++) {
        for (int q=0; q < numOrb; q++) {
            int pq = compoundIndex(p, q);
            for (int r=0; r < numOrb; r++) {
                for (int s=0; s < numOrb; s++) {
                    int rs = compoundIndex(r, s);
                    int pqrs = compoundIndex(pq,rs);
                    TEI_MOspin[p][q][r][s] = (*TEI_MO)(pqrs) * (p%2 == r%2) * (q%2 == s%2);
                }
            }
        }
    }
}

/**
 AO to MO integral transformation using a single N^8 step
 Both Two Electron Integrals are stored in arrays, taking advantage of
 Permuational Symmetry

 runtime : O(N^8)
 **/
void atomicToMolecularN8(MatrixXd *TEI_MO, MatrixXd *TEI_AO, MatrixXd *C0) {

    int numOrb = (*C0).rows();
    int i, j, k, l, ijkl;
    int p, q, r, s, pqrs;

    setzero(TEI_MO);

    for(i=0, ijkl=0; i < numOrb; i++) {
    for(j=0; j <= i; j++) {
    for(k=0; k <= i; k++) {
    for(l=0; l <= (i==k ? j : k); l++, ijkl++) {

    for(p=0; p < numOrb; p++) { // Over all orbitals
    for(q=0; q < numOrb; q++) { // Over all orbitals

        for(r=0; r < numOrb; r++) { // Over all orbitals
        for(s=0; s < numOrb; s++) { // Over all orbitals

        pqrs = compoundIndex(compoundIndex(p,q),compoundIndex(r,s));

        (*TEI_MO)(ijkl) += (*C0)(p,i) * (*C0)(q,j) * (*C0)(r,k) * (*C0)(s,l) * (*TEI_AO)(pqrs);

        }}
    }}}}}}
  }

/**
 AO to MO integral transformation using a single N^8 step
 Both Two Electron Integrals are stored in arrays, taking advantage of
 Permuational Symmetry

 runtime : O(N^5)
 **/
void atomicToMolecularN5(MatrixXd *TEI_MO, MatrixXd *TEI_AO, MatrixXd *C0) {

      if (C0 == NULL || TEI_AO == NULL || TEI_MO == NULL){
        throw invalid_argument( "received null pointer in atomicToMolecularN5" );
      }

      int numOrb = (*C0).rows();
      MatrixXd* X = new MatrixXd(numOrb, numOrb);
      setzero(X);

      MatrixXd* Temp = new MatrixXd(NUM_MUL, NUM_MUL);
      setzero(Temp);

      int i, j, ij, k, l, kl;

      for (i = 0, ij = 0; i < numOrb; i++) {
          for (j = 0; j <= i; j++, ij++) {
              for (k = 0, kl = 0; k < numOrb; k++) {
                  for (l = 0; l <= k; l++, kl++) {
                      (*X)(k,l) = (*X)(l,k) = (*TEI_AO)(compoundIndex(ij,kl));
                    }
              }
              (*X) = (*C0).transpose() * (*X) * (*C0);
              for (k = 0, kl = 0; k < numOrb; k++) {
                  for (l = 0; l <= k; l++, kl++) {
                      (*Temp)(kl,ij) = (*X)(k,l);
                  }
              }
          }
      }

      setzero(TEI_MO);

      for (k = 0, kl = 0; k < numOrb; k++) {
          for (l = 0; l <= k; l++, kl++) {
              for (i = 0, ij = 0; i < numOrb; i++) {
                  for (j = 0; j <= i; j++, ij++) {
                      (*X)(i,j) = (*X)(j,i) = (*Temp)(kl,ij);
                  }
              }
              (*X) = (*C0).transpose() * (*X) * (*C0);
              for (i = 0, ij = 0; i < numOrb; i++){
                  for (j = 0; j <= i; j++, ij++){
                      (*TEI_MO)(compoundIndex(kl,ij)) = (*X)(i,j);
                  }
             }
          }
      }

      delete Temp;
      delete X;

  }

/**
  Creates the Spin Orbital Fock Matrix, given a set of TEI
  in the Spin Orbital Basis and the one elctron hamiltonian

  @NOTE:
    requires all occupied orbitals before virtual orbitals

  runtime : O(N_occ*N^2)
  **/
void spinOrbitalFock(MatrixXd *FSO, double ****TEI_MOspin, MatrixXd *h, int numOcc) {

    // Fij = hij + Sum over occupied orbitals <pm||qm>
    // where <pm||qm> = <pm|qm> - <pm|mq> = [pq|mm] - [pm|mq]
    if (FSO == NULL || TEI_MOspin == NULL || h == NULL){
        throw invalid_argument( "received null pointer in spinOrbitalFock" );
    }

    int numOrb = (*h).rows();

    for (int p=0; p < numOrb; p++) {
      for (int q=0; q < numOrb; q++) {

        (*FSO)(p,q) = (*h)(p,q); // one electron hamiltonians

        for(int m=0; m < numOcc; m++) {
          (*FSO)(p,q) += getExchangeIntegral(TEI_MOspin, p, m, q, m);
        }
      }
    }
  }

double MP2_Energy(MatrixXd *TEI_MO, MatrixXd *E, int nElectrons) {
    // E are orbital Energy values

    if (E == NULL || TEI_MO == NULL){
        throw invalid_argument( "received null pointer in MP2_Energy" );
    }

    double EMP2 = 0;
    int numOrb = (*E).rows();
    int nOcc = nElectrons/2; // closed shell
    int ia, ja, jb, ib, iajb, ibja;

    for (int i = 0; i < nOcc; i++){
        for (int a = nOcc; a < numOrb; a++){
            ia = compoundIndex(i,a);
            for (int j = 0; j < nOcc ; j++){
                ja = compoundIndex(j,a);
                for (int b = nOcc; b < numOrb; b++){
                    jb = compoundIndex(j,b);
                    ib = compoundIndex(i,b);
                    iajb = compoundIndex(ia,jb);
                    ibja = compoundIndex(ib,ja);
                    EMP2 += (*TEI_MO)(iajb) * ( 2.0 * (*TEI_MO)(iajb) - (*TEI_MO)(ibja) ) / ((*E)(i) + (*E)(j) - (*E)(a) - (*E)(b));
                }
            }
        }
    }
    return EMP2;
}

void copyMatrix(MatrixXd *M, MatrixXd *N) {

    if (M == NULL || N == NULL){
        throw invalid_argument( "received null pointer in copyMatrix" );
    }

    if ((*M).rows() != (*N).rows() || (*M).cols() != (*N).cols()) {
        string msg = "Copy matrix must be of the same dimension : \n";
        char stemp[100] = "";
        char stemp2[100] = "";
        snprintf(stemp, 100, "first of dim (%ld, %ld) \n", (*M).rows() , (*M).cols());
        snprintf(stemp2, 100, "second of dim (%ld, %ld)", (*N).rows() , (*N).cols());
        cout << msg + stemp + stemp2 << endl;
        exit(-1);
    }

    for (int i = 0; i < (*M).rows(); i++){
        for (int j = 0; j < (*M).cols(); j++){
            (*N)(i,j) = (*M)(i,j);
        }
    }
}

double getExchangeIntegral(MatrixXd *TEI, int i, int j, int a, int b) {

    if (TEI == NULL){
        throw invalid_argument( "received null pointer in getExchangeIntegral" );
    }

    int iajb = compoundIndex(compoundIndex(i,a),compoundIndex(j,b));
    int ibja = compoundIndex(compoundIndex(i,b),compoundIndex(j,a));
    return (*TEI)(iajb) - (*TEI)(ibja);
}

double getExchangeIntegral(double ****TEI, int i, int j, int a, int b) {

    if (TEI == NULL){
        throw invalid_argument( "received null pointer in getExchangeIntegral" );
    }
    return TEI[i][a][j][b] - TEI[i][b][j][a];
}

void printMatrix(MatrixXd *M, string name){

    if (M == NULL){
        throw invalid_argument( "attemping to print an unset matrix" );
    }

    cout << "--------------------------------------------------------------------------------" << endl;
    cout << name << " : " << (*M) << endl;
}