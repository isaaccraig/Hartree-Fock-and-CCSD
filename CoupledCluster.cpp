//
// Created by Isabel Craig on 3/26/18.
//
// TODO :
// Necessary contractions and symmetry treatment as described in
// reference 1 : J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett,
// J. Chem. Phys. volume 94, pp. 4334-4345 (1991).
//
//

#include "CoupledCluster.hpp"

using namespace std;
using namespace Eigen;

#define OCCUPIED_RANGE(i) (int i = 0; i < numOcc; i++)
#define UNOCCUPIED_RANGE(a) (int a = numOcc; a < numOrb; a++)
#define ij (i * numOrb + j)
#define ab (a * numOrb + b)

CCSD::CCSD(HartreeFock HF, double tolerance_e) {

    this->tol_e = tolerance_e;
    this->numOrb = HF.numBasisFunc;
    this->numOcc = 0.5 * HF.numElectrons;
    this->numUnocc = numOrb - numOcc;

    this->tau = new MatrixXd(numOrb * numOrb, numOrb * numOrb);
    this->tautilda = new MatrixXd(numOrb * numOrb, numOrb * numOrb);
    this->F = new MatrixXd(numOrb, numOrb);
    this->Fock = new MatrixXd(numOrb, numOrb);
    this->TEI = allocate4DMatrix(numOrb);
    this->W = allocate4DMatrix(numOrb);

    molecularToMolecularSpin(this->TEI, HF.TEI_MO, numOrb);
    spinOrbitalFock(this->Fock, this->TEI, HF.Hcore, numOcc);

    setDenominatorArrays();
    setInitialAmplitudes(HF.orbitalE);
    setExcitation();
    setIntermediates();
    this->E = correlationEnergy();
    this->prevE = 100;

}

void CCSD::printState(){

    printMatrix(tau, "tau");
    printMatrix(tautilda, "tautilda");
    printMatrix(tau, "tau");
    printMatrix(F, "F");
    printMatrix(Fock, "Fock");
    printMatrix(D_1d, "D1d");
    printMatrix(D_2d, "D2d");

}

CCSD::~CCSD(){
    delete tau;
    delete tautilda;
    delete F;
    delete Fock;
    delete4DMatrix(TEI, numOrb);
    delete4DMatrix(W, numOrb);
}

void CCSD::Iterate(){
    cout << "Iteration : \t\t" << "Energy : \t\t" << endl;
    int it = 0;
    while (!eConverg()){
        printState();
        setExcitation();
        setIntermediates();
        updateAmplitudes();
        prevE = E;
         E = correlationEnergy();
        cout << E << "\t\t" << it << endl;
    }
}

void CCSD::setExcitation(){
    // effective doubles
    // see equations 9, 10 in reference 1

    copyMatrix(t_double, tau);
    copyMatrix(t_double, tautilda);

    for OCCUPIED_RANGE(i) {
        for UNOCCUPIED_RANGE(a) {
            for OCCUPIED_RANGE(j) {
                for UNOCCUPIED_RANGE(b) {
                    double val = (*t_single)(i,a) * (*t_single)(j,b) - (*t_single)(i,b) * (*t_single)(j,a);
                    (*tautilda)(numOrb * i + j, numOrb * a + b) += 0.5 * val;
                    (*tau)(numOrb * i + j, numOrb * a + b) += val;
                }
            }
        }
    }

}

void CCSD::setIntermediates(){

    // must first set excitation operators tau and tau tilda
    // using CCSD::setExcitation

    // sets F (one index intermediate)

    for UNOCCUPIED_RANGE(a) {
        for UNOCCUPIED_RANGE(e) {
            // p unoccupied orbital
            // Equation 3 in Reference 1
            (*F)(a, e) = (1 - (a == e)) * (*Fock)(a, e);
            for OCCUPIED_RANGE(m) {
                (*F)(a, e) -= 0.5 * (*t_single)(m, a) * (*Fock)(m, e);
                for UNOCCUPIED_RANGE(f) {
                    (*F)(a, e) += (*t_single)(m, f) * getExchangeIntegral(TEI, m, a, f, e);
                    for (int n = 0; n < numOrb; n++) {
                        (*F)(a, e) -= 0.5 * getTauTilda(m, n, a, f) * getExchangeIntegral(TEI, m, n, e, f);
                    }
                }
            }

        }
    }

    for OCCUPIED_RANGE(m) {
        for OCCUPIED_RANGE(i) {
            // q occupied orbital
            // Equation 4 in Reference 1
            (*F)(m, i) = (1 - (m == i)) * (*Fock)(m, i);
            for UNOCCUPIED_RANGE(e) {
                (*F)(m, i) += 0.5 * (*t_single)(i, e) * (*Fock)(m, e);
                for OCCUPIED_RANGE(n) {
                    (*F)(m, i) += (*t_single)(n, e) * getExchangeIntegral(TEI, m, n, i, e);
                    for UNOCCUPIED_RANGE(f) {
                        (*F)(m, i) += 0.5 * getTauTilda(i, n, e, f) * getExchangeIntegral(TEI, m, n, e, f);
                    }
                }
            }
        }
    }

    for OCCUPIED_RANGE(m) {
        for UNOCCUPIED_RANGE(e) {
            // Equation 5 in Reference 1
            (*F)(m, e) = (*Fock)(m, e);
            for OCCUPIED_RANGE(n) {
                for UNOCCUPIED_RANGE(f) {
                    (*F)(m, e) += (*t_single)(n, f) * getExchangeIntegral(TEI, m, n, e, f);
                }
            }
        }
    }

    for OCCUPIED_RANGE(m) {
        for OCCUPIED_RANGE(n) {
            for OCCUPIED_RANGE(i) {
                for OCCUPIED_RANGE(j) {
                    // Equation 6 in Reference 1
                    W[m][n][i][j] = getExchangeIntegral(TEI, m, n, i, j);
                    for UNOCCUPIED_RANGE(e) {
                        W[m][n][i][j] +=  (*t_single)(j, e) * getExchangeIntegral(TEI, m, n, i, e);
                        W[m][n][i][j] -=  (*t_single)(i, e) * getExchangeIntegral(TEI, m, n, j, e);
                        for UNOCCUPIED_RANGE(f) {
                            W[m][n][i][j] += 0.25 * getTau(i, j, e, f) * getExchangeIntegral(TEI, m, n, e, f);
                        }
                    }
                }
            }
        }
    }

    for UNOCCUPIED_RANGE(a) {
        for UNOCCUPIED_RANGE(b) {
            for UNOCCUPIED_RANGE(e) {
                for UNOCCUPIED_RANGE(f) {
                    // Equation 7 in Reference 1
                    W[a][b][e][f] = getExchangeIntegral(TEI, a, b, e, f);
                    for OCCUPIED_RANGE(m) {
                        W[a][b][e][f] -=  (*t_single)(m, b) * getExchangeIntegral(TEI, a, m, e, f);
                        W[a][b][e][f] +=  (*t_single)(m, a) * getExchangeIntegral(TEI, b, m, e, f);
                        for OCCUPIED_RANGE(n) {
                            W[a][b][e][f] += 0.25 * getTau(m, n, a, b) * getExchangeIntegral(TEI, m, n, e, f);
                        }
                    }
                }
            }
        }
    }

    for OCCUPIED_RANGE(m) {
        for UNOCCUPIED_RANGE(b) {
            for UNOCCUPIED_RANGE(e) {
                for UNOCCUPIED_RANGE(j) {
                    // Equation 8 in Reference 1
                        W[m][b][e][j] = getExchangeIntegral(TEI, m, b, e, j);
                    for UNOCCUPIED_RANGE(f) {
                        W[m][b][e][j] += (*t_single)(j, f) * getExchangeIntegral(TEI, m, b, e, f);
                    }
                    for OCCUPIED_RANGE(n) {
                        W[m][b][e][j] -= (*t_single)(n, b) * getExchangeIntegral(TEI, m, n, e, j);
                    }
                    for OCCUPIED_RANGE(n) {
                        for UNOCCUPIED_RANGE(f) {
                            W[m][b][e][j] -= getExchangeIntegral(TEI, m, n, e, f) *
                                    (0.5 * getDoubleAmplitude(j, n, f, b) - (*t_single)(j, f) * (*t_single)(n, b));
                        }
                    }
                }
            }
        }
    }
}

void CCSD::setDenominatorArrays(){
    // see equations 12, 13 in reference 1
    D_1d = new MatrixXd(numOrb, numOrb);
    D_2d = new MatrixXd(numOrb * numOrb, numOrb * numOrb);

    for OCCUPIED_RANGE(i) {
        for UNOCCUPIED_RANGE(a) {
            (*D_1d)(i,a) = (*Fock)(i,i) - (*Fock)(a,a);
        }
    }
    for OCCUPIED_RANGE(i) {
        for UNOCCUPIED_RANGE(a) {
            for OCCUPIED_RANGE(j) {
                for UNOCCUPIED_RANGE(b) {
                    (*D_2d)(i*numOrb + j, a*numOrb + b) = (*Fock)(i, i) + (*Fock)(j, j) - (*Fock)(a, a) - (*Fock)(b, b);
                }
            }
        }
    }
}

void CCSD::updateAmplitudes(){
    // see equations 1, 2 in reference 1

    auto *M_1d = new MatrixXd(numOrb, numOrb);

    for OCCUPIED_RANGE(i) {
        for UNOCCUPIED_RANGE(a) {
            (*M_1d)(i,a) = (*Fock)(i, a);
            for UNOCCUPIED_RANGE(e) {
                (*M_1d)(i,a) += (*t_single)(e, i) * (*F)(a,e);
            }
            for OCCUPIED_RANGE(m) {
                (*M_1d)(i,a) -= (*t_single)(m, a) * (*F)(m,i);
                for UNOCCUPIED_RANGE(e) {
                    (*M_1d)(i,a) += getDoubleAmplitude(i, m, a, e) * (*F)(m,e);
                    for UNOCCUPIED_RANGE(f) {
                        (*M_1d)(i,a) -= 0.5 * getDoubleAmplitude(i, m, e, f) * getExchangeIntegral(TEI, m, a, e, f);
                    }
                    for OCCUPIED_RANGE(n) {
                        (*M_1d)(i,a) -= 0.5 * getDoubleAmplitude(a, e, m, n) * getExchangeIntegral(TEI, m, n, e, i);
                    }
                }
            }
            for UNOCCUPIED_RANGE(f) {
                for OCCUPIED_RANGE(n) {
                    (*M_1d)(i,a) -= (*t_single)(n, f) * getExchangeIntegral(TEI, n, a, i, f);
                }
            }
        }
    }

    auto *M_2d = new MatrixXd(numOrb*numOrb, numOrb*numOrb);

    for OCCUPIED_RANGE(i) {
        for OCCUPIED_RANGE(j) {
            for UNOCCUPIED_RANGE(a) {
                for UNOCCUPIED_RANGE(b) {

                    (*M_2d)(ij, ab) = getExchangeIntegral(TEI, i, j, a, b);

                    for UNOCCUPIED_RANGE(e) {
                        for UNOCCUPIED_RANGE(f) {
                            (*M_2d)(ij, ab) += 0.5 * getTau(i, j, e, f) * getW(a, b, e, f);
                        }
                    }
                    for OCCUPIED_RANGE(m) {
                        for OCCUPIED_RANGE(n) {
                            (*M_2d)(ij, ab) += 0.5 * getTau(m, n, a, b) * getW(m, n, i, j);
                        }
                    }
                    for UNOCCUPIED_RANGE(e) {

                        (*M_2d)(ij, ab) += (*t_single)(i, e) * getExchangeIntegral(TEI, a, b, e, j);
                        (*M_2d)(ij, ab) -= (*t_single)(j, e) * getExchangeIntegral(TEI, a, b, e, i);

                        int val = 0;
                        for OCCUPIED_RANGE(m) {
                            val += (*t_single)(m, b) * (*F)(m,e);
                        }

                        (*M_2d)(ij, ab) += getDoubleAmplitude(i, j, a, e) * ( (*F)(b,e) - 0.5 * val);
                        (*M_2d)(ij, ab) -= getDoubleAmplitude(i, j, b, e) * ( (*F)(a,e) - 0.5 * val);

                    }
                    for OCCUPIED_RANGE(m) {

                        (*M_2d)(ij, ab) -= (*t_single)(m, a) * getExchangeIntegral(TEI, m, b, i, j);
                        (*M_2d)(ij, ab) += (*t_single)(m, b) * getExchangeIntegral(TEI, m, a, i, j);

                        int val = 0;
                        for UNOCCUPIED_RANGE(e) {
                            val += (*t_single)(j, e) * (*F)(m,e);
                        }

                        (*M_2d)(ij, ab) -= getDoubleAmplitude(i, m, a, b) * ( (*F)(m,j) + 0.5 * val);
                        (*M_2d)(ij, ab) += getDoubleAmplitude(j, m, a, b) * ( (*F)(m,i) + 0.5 * val);

                    }

                    for UNOCCUPIED_RANGE(e) {
                        for OCCUPIED_RANGE(m) {

                            (*M_2d)(ij, ab) += getDoubleAmplitude(i, m, a, e) * getW(m, b, e, j);
                            (*M_2d)(ij, ab) += getDoubleAmplitude(j, m, b, e) * getW(m, a, e, i);
                            (*M_2d)(ij, ab) -= getDoubleAmplitude(i, m, b, e) * getW(m, a, e, j);
                            (*M_2d)(ij, ab) -= getDoubleAmplitude(j, m, a, e) * getW(m, b, e, i);

                            (*M_2d)(ij, ab) -= (*t_single)(i, e) * (*t_single)(m, a) * getExchangeIntegral(TEI, m, b, e, j);
                            (*M_2d)(ij, ab) -= (*t_single)(j, e) * (*t_single)(m, b) * getExchangeIntegral(TEI, m, a, e, i);
                            (*M_2d)(ij, ab) += (*t_single)(j, e) * (*t_single)(m, a) * getExchangeIntegral(TEI, m, b, e, i);
                            (*M_2d)(ij, ab) += (*t_single)(i, e) * (*t_single)(m, b) * getExchangeIntegral(TEI, m, a, e, j);
                        }
                    }
                }
            }
        }
    }

    (*t_single) = (*M_1d) * (*D_1d).inverse();
    (*t_double) = (*M_2d) * (*D_2d).inverse();

    delete M_1d;
    delete M_2d;

}

bool CCSD::eConverg(){
    // checks for convergence of the energy value
    double delE = (prevE - E);
    return (delE < tol_e);
}

double CCSD::correlationEnergy(){
    double E = 0;

    for OCCUPIED_RANGE(i){
        for UNOCCUPIED_RANGE(a) {
            E += (*t_single)(i, a) * (*Fock)(i, a);
        }
    }

    for OCCUPIED_RANGE(i){
        for UNOCCUPIED_RANGE(a){
            for OCCUPIED_RANGE(j){
                for UNOCCUPIED_RANGE(b){
                    E += 0.25 * getDoubleAmplitude(i,j,a,b) * getExchangeIntegral(TEI, i, j, a, b);
                }
            }
        }
    }

    for OCCUPIED_RANGE(i){
        for UNOCCUPIED_RANGE(a){
            for OCCUPIED_RANGE(j){
                for UNOCCUPIED_RANGE(b){
                    E += 0.5 * (*t_single)(i, a) * (*t_single)(j, b) * getExchangeIntegral(TEI, i, j, a, b);
                }
            }
        }
    }

    return E;

}

void CCSD::setInitialAmplitudes(MatrixXd *energies){

    t_single = new MatrixXd(numOrb, numOrb);
    t_double = new MatrixXd(numOrb * numOrb, numOrb * numOrb);

    for OCCUPIED_RANGE(i){
        for UNOCCUPIED_RANGE(a){
            (*t_single)(i, a) = 0;
            for OCCUPIED_RANGE(j){
                for UNOCCUPIED_RANGE(b){
                    (*t_double)(numOrb * i + j, numOrb * a + b) = getExchangeIntegral(TEI, i, j, a, b) /
                            ((*energies)(i) + (*energies)(i) - (*energies)(a) - (*energies)(b));
                }
            }
        }
    }

}

double CCSD::getDoubleAmplitude(int i, int j, int a, int b){
    return (*t_double)( numOrb * i + j, numOrb * a + b);
}

double CCSD::getTauTilda(int i, int j, int a, int b){
    return (*tautilda)( numOrb * i + j, numOrb * a + b);
}

double CCSD::getTau(int i, int j, int a, int b){
    return (*tau)( numOrb * i + j, numOrb * a + b);
}

double CCSD::getW(int i, int j, int a, int b){
    return W[i][j][a][b];
}

double CCSD::MP2Energy(){
    double emp2 = 0;
    for OCCUPIED_RANGE(i){
        for UNOCCUPIED_RANGE(a){
            for OCCUPIED_RANGE(j){
                for UNOCCUPIED_RANGE(b){
                    emp2 += getDoubleAmplitude(i,j,a,b) * getExchangeIntegral(TEI, i, j, a, b);
                }
            }
        }
    }
    cout << "CCSD Estimated EMP2 :" << emp2 << endl;
    return emp2;
}