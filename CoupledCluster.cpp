//
// Created by Isabel Craig on 3/26/18.
//

#include "CoupledCluster.hpp"
#include "Read.hpp"
#include "QuantumUtils.hpp"
using namespace std;
using namespace Eigen;

CoupledCluster::CoupledCluster(MatrixXd *h, MatrixXd *tei){
    molecularToMolecularSpin(this->TEI, tei);
    SpinOrbitalFock(this->Fock, this->TEI, h);
}

