

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
  int i, j, k, l, ij, kl, ji, lk;
  double val;

  while(fscanf(input, "%d %d %d %d %lf", &i, &j, &k, &l, &val) != EOF) {

    ij = INDEX(i-1,j-1);
    kl = INDEX(k-1,l-1);
    ji = INDEX(j-1,i-1);
    lk = INDEX(l-1,k-1);
      
 /**
    if (ijkl > (*TEI).size()) {
        cout << "Error : Matrix Size Exceded" << endl;
        cout << (*TEI).size() << " = size" << endl;
        cout << ijkl << " = index" << endl;
        exit(-1);
    } **/
      
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

void READIN::Test_Mulliken(MullikenMatrix *TEI) {
    
    int ij, ji, kl, lk;
    bool val;
    
    for (int i = 0; i < NUM_ORB ; i++) {
        for (int j = 0; j < NUM_ORB ; j++) {
            for (int k = 0; k < NUM_ORB ; k++) {
                for (int l = 0; l < NUM_ORB ; l++) {
                    
                    ij = INDEX(i,j);
                    ji = INDEX(j,i);
                    kl = INDEX(k,l);
                    lk = INDEX(l,k);
                    
                    val = ((*TEI)(INDEX(ij,kl)) ==
                           (*TEI)(INDEX(ji,kl)) ==
                           (*TEI)(INDEX(ij,lk)) ==
                           (*TEI)(INDEX(ji,lk)) ==
                           (*TEI)(INDEX(kl,ij)) ==
                           (*TEI)(INDEX(kl,ji)) ==
                           (*TEI)(INDEX(lk,ji)) ==
                           (*TEI)(INDEX(lk,ij)));
                    
                    if (!val) {
                        cout << "MULIKEN FAILED" << endl;
                        cout << ij << " " << ji << " " << kl << " " << lk << " " <<endl;
                        exit(-1);
                    }
                }
            }
        }
    }
    
}

