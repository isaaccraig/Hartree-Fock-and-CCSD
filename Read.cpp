

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

void READIN::Mulliken(const char *filename, double (*TEI)[NUM_ORB][NUM_ORB][NUM_ORB][NUM_ORB]){

  FILE *input;
  input = fopen(filename, "r");
  int i, j, k, l;
  double val;

  for (int i = 0; i < NUM_ORB; i++) {
      for (int j = 0; j < NUM_ORB; j++) {
          for (int k = 0; k < NUM_ORB; k++) {
              for (int l = 0; l < NUM_ORB; l++) {
                  (*TEI)[i][j][k][l] = FILLER;
              }
          }
      }
  }


  while(fscanf(input, "%d %d %d %d %lf", &i, &j, &k, &l, &val) != EOF) {

    i--;
    j--;
    k--;
    l--;

    (*TEI)[i][j][k][l] = val;
    (*TEI)[j][i][k][l] = val;
    (*TEI)[i][j][l][k] = val;
    (*TEI)[j][i][l][k] = val;

    (*TEI)[k][l][i][j] = val;
    (*TEI)[k][l][j][i] = val;
    (*TEI)[l][k][i][j] = val;
    (*TEI)[l][k][j][i] = val;

  }
  fclose(input);
}

void READIN::Test_Mulliken(double (*TEI)[NUM_ORB][NUM_ORB][NUM_ORB][NUM_ORB]) {

    bool val;

    for (int i = 0; i < NUM_ORB ; i++) {
        for (int j = 0; j < NUM_ORB ; j++) {
            for (int k = 0; k < NUM_ORB ; k++) {
                for (int l = 0; l < NUM_ORB ; l++) {

                  val = ( (*TEI)[i][j][k][l] ==
                          (*TEI)[j][i][k][l] ==
                          (*TEI)[i][j][l][k] ==
                          (*TEI)[j][i][l][k] ==
                          (*TEI)[k][l][i][j] ==
                          (*TEI)[k][l][j][i] ==
                          (*TEI)[l][k][i][j] ==
                          (*TEI)[l][k][j][i]);

                    if (!val) {
                        cout << "MULIKEN FAILED" << endl;
                        exit(-1);
                    }
                }
            }
        }
    }
}
