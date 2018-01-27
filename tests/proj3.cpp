#include "proj3.hpp"

using namespace std;
using namespace Eigen;

//Read in matrix from file.
void proj3::readMatrix2d(const char* fname, MatrixXd* m)
{
	FILE *input;
 	input = fopen(fname, "r");
 	int i,j;
 	double x;
 	while( !feof(input) )
 	{
 		fscanf(input, "%d %d %lf", &i, &j, &x);
    	(*m)(i-1, j-1) = x;
 	}
}

//Reflect matrix acros diagonal assuming it is square - can choose upper triangular half to reflect with
//of lower diagonal half: 0-upper, 1-lower
void proj3::reflectMatrix2d(MatrixXd* m, int mode)
{
	int size = m->rows();
	if(mode)
	{
		//Traverse down column
		for(int i = 0; i < size; i++)
		{
			//Traverse across row
			for(int j = i + 1; j < size; j++)
			{
				(*m)(i, j) = (*m)(j, i);
			}
		}
	}
	else
	{
		//Traverse across row
		for(int i = 0; i < size; i++)
		{
			//Traverse down column
			for(int j = i + 1; j < size; j++)
			{
				(*m)(j, i) = (*m)(i, j);
			}
		}
	}
}

//Create 4D Matrix
MatrixXd*** proj3::create4DMatrix(MatrixXd*** m, int sizes[])
{
	m = new MatrixXd** [sizes[0]];
	for (int i = 0; i < sizes[0]; i++) 
	{
        m[i] = new MatrixXd* [sizes[1]];
        for (int j = 0; j < sizes[1]; j++) 
            m[i][j] = new (MatrixXd)(MatrixXd::Ones(sizes[2],sizes[3]).array() * 0.0);
    }
    return m;
}


//Read in 4D matrix in Mulliken notation. I am also assuming it is symmetric.
//MAKE SURE THE DIMENSIONS ARE EQUAL IN EVERY AXIS OR YOU WILL SEGFAULT.
void proj3::readMatrix4d(const char* fname, MatrixXd*** m)
{
	FILE *input;
 	input = fopen(fname, "r");
 	int i,j,k,l;
 	double x;
 	while( !feof(input) )
 	{
 		fscanf(input, "%d %d %d %d %lf", &i, &j, &k, &l, &x);
    	(*m[i-1][j-1])(k-1, l-1) = x;	//ijkl
    	(*m[j-1][i-1])(k-1, l-1) = x;	//jikl
    	(*m[i-1][j-1])(l-1, k-1) = x;	//ijlk
    	(*m[j-1][i-1])(l-1, k-1) = x;	//jilk
    	(*m[k-1][l-1])(i-1, j-1) = x;	//klij
    	(*m[l-1][k-1])(i-1, j-1) = x;	//lkij
    	(*m[k-1][l-1])(j-1, i-1) = x;	//klji
    	(*m[l-1][k-1])(j-1, i-1) = x;	//lkji
 	}
} 

//Print an error message and show faulty matrix comparison.
void proj3::printErr(char* msg, Eigen::MatrixXd src, Eigen::MatrixXd bad)
{
	cout << msg << endl << endl;
	cout << "Actual answer:" << endl << src << endl << endl;
	cout << "Your answer:" << endl << bad << endl << endl;
	assert(false);
}