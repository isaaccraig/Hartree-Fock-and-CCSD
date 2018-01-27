#include <iostream>
#include <string>
#include <list>
#include <assert.h>

#include "Eigen/Dense"
#include "Eigen/Core"

#include "../HartreeFock.hpp"
#include "../Read.hpp"

//Basic helper constructs
class proj3
{
public:
	int power(int, int);

	//2D Matrix operations
	void readMatrix2d(const char*, Eigen::MatrixXd*);
	void reflectMatrix2d(Eigen::MatrixXd*, int);

	//4D Matrices
	Eigen::MatrixXd*** create4DMatrix(Eigen::MatrixXd***, int[]);
	void readMatrix4d(const char*, Eigen::MatrixXd***);

	//Print Error in comparing matrices
	void printErr(char* msg, Eigen::MatrixXd, Eigen::MatrixXd);
};