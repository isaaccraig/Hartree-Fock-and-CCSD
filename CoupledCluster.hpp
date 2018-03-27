//
// Created by Isabel Craig on 3/26/18.
//

#ifndef HARTREEFOCKMP2_COUPLEDCLUSTER_H
#define HARTREEFOCKMP2_COUPLEDCLUSTER_H

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

class CoupledCluster{

private:
    Eigen::MatrixXd *Fock;
    Eigen::MatrixXd *TEI;

public:
    CoupledCluster();
    ~CoupledCluster();

};











#endif //HARTREEFOCKMP2_COUPLEDCLUSTER_H
