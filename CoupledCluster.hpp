//
// Created by Isabel Craig on 3/26/18.
//

#ifndef HARTREEFOCKMP2_COUPLEDCLUSTER_H
#define HARTREEFOCKMP2_COUPLEDCLUSTER_H

#include "Read.hpp"
#include "QuantumUtils.hpp"
#include "HartreeFock.hpp"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

class CCSD{

private:
    Eigen::MatrixXd *Fock;
    Eigen::MatrixXd *TEI;

    Eigen::MatrixXd *tau;
    Eigen::MatrixXd *tautilda;
    Eigen::MatrixXd *t_double;
    Eigen::MatrixXd *t_single;
    Eigen::MatrixXd *D_1d;
    Eigen::MatrixXd *D_2d;
    double ****W;
    Eigen::MatrixXd *F;

    double tol_e;
    double prevE;
    double E;

    int numOcc;
    int numOrb;
    int numUnocc;

    void setExcitation();
    void setIntermediates();
    void setDenominatorArrays();
    bool eConverg();
    void updateAmplitudes();
    void setInitialAmplitudes(Eigen::MatrixXd *energies);
    double getDoubleAmplitude(int i, int j, int a, int b);
    double getTauTilda(int i, int j, int a, int b);
    double getTau(int i, int j, int a, int b);
    double getW(int i, int j, int a, int b);

public:
    CCSD(HartreeFock HF, double tol_e);
    double MP2Energy();
    void printState();
    double correlationEnergy();
    void Iterate();
    ~CCSD();

};

#endif //HARTREEFOCKMP2_COUPLEDCLUSTER_H
