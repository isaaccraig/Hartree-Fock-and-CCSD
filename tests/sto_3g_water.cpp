#include "proj3.hpp"

using namespace std;
using namespace Eigen;

/* -- Test suite for sto_3g_water, to make tests for the others, you can actually 
just use the same code with different .dat files. -- */

class sto_3g_water : proj3
{
public:
    sto_3g_water(double, double);
    void destroy(void);

    void test1(void);
    void test2(void);
    void test3(void);
    void test4(void);
    void test5(void);
private:
    HartreeFock* hfTestObj;
};

/* --- Test Suite tests  --- */
//Test 1 - Check nuclear repulsion energy output
void sto_3g_water::test1(void)
{
	cout << "Runing Test 1." << endl;
	cout << "Running test1: Checking nuclear repulsion energy..." << endl;
	assert(abs(8.0023670618 - hfTestObj->enuc) < 1e-7);
	cout << "Test 1 passed." << endl << endl;
}

//Test 2 - Checking Integral matrices
void sto_3g_water::test2(void)
{
	cout << "Running Test 2." << endl;
	cout << "Running test2: Checking overlap matrix..." << endl;
	MatrixXd test(7, 7);
	readMatrix2d("data/overlap.dat", &test);
	reflectMatrix2d(&test, 1);
	if(!test.isApprox( hfTestObj->S ))
		printErr( (char*) "Overlap matrices are not correct...", test, hfTestObj->S );
	
	cout << "Running test2: Checking kinetic energy matrix..." << endl;
	readMatrix2d("data/kinetic.dat", &test);
	reflectMatrix2d(&test, 1);
	if(!test.isApprox( hfTestObj->T ))
		printErr( (char*) "KE matrix is not correct...", test, hfTestObj->T );

	cout << "Running test2: Checking nuclear attraction matrix..." << endl;
	readMatrix2d("data/anuc.dat", &test);
	reflectMatrix2d(&test, 1);
	if(!test.isApprox( hfTestObj->V ))
		printErr( (char*) "Nuclear attraction matrix is not correct...", test, hfTestObj->V );

	cout << "Test 2 passed." << endl << endl;
}

//Test 3 - Checking Two electron integrals
void sto_3g_water::test3(void)
{
	cout << "Running Test 3." << endl;
	cout << "Running test3: Checking Electron Repulsion 4D Matrix..." << endl;
	MatrixXd*** test;
	int size = 7;
	int sizes[4] = {size, size, size, size};
	test = create4DMatrix(test, sizes);
	readMatrix4d("data/eri.dat", test);

	MullikenMatrix TEI = hfTestObj->TEI;

	for( int i = 0; i < size; i++ )
	{
		for( int j = 0; j < size; j++ )
		{
			for( int k = 0; k < size; k++ )
			{
				for( int l = 0; l < size; l++ )
				{
					int ij = INDEX(i,j);
				    int kl = INDEX(k,l);
				    int ji = INDEX(j,i);
				    int lk = INDEX(l,k);
					if( abs( (*test[i][j])(k, l)  - TEI(INDEX(ij,kl)) )  > 1e-4 )
					{
						cout << "Your 4D Matrix is wrong at index: i = " << i << ", j = " << j << ", k = " << k << ",l = " << l << endl;
						assert(false);
					}
				}
			}
		}
	}

	cout << "Test 3 passed." << endl << endl;
}

//Test 4 - Testing Orthogonal Matrix
void sto_3g_water::test4(void)
{
	cout << "Running Test 4." << endl;
	cout << "Running test4: Checking Orthogonality Matrix..." << endl;
	MatrixXd test(7, 7);
	readMatrix2d("data/orthog.dat", &test);
	reflectMatrix2d(&test, 1);
	if(!test.isApprox( hfTestObj->SOM, 1e-5 ))
		printErr( (char*) "Orthogonality matrix is not correct...", test, hfTestObj->SOM );
	cout << "Test 4 passed." << endl << endl;
}

//Test 5 - Testing Initial Density Matrix
void sto_3g_water::test5(void)
{
	cout << "Running Test 5." << endl;
	cout << "Running test4: Checking Initial Density Matrix..." << endl;
	MatrixXd test(7, 7);
	readMatrix2d("data/fock_t.dat", &test);
	reflectMatrix2d(&test, 1);
	if(!test.isApprox( hfTestObj->F0, 1e-5 ))
		printErr( (char*) "Inital fock matrix is not correct...", test, hfTestObj->F0 );
	cout << "Test 5 passed." << endl << endl;

	cout << "Running test4: Checking Initial MO Matrix..." << endl;
	readMatrix2d("data/mo_i.dat", &test);
	if(!test.isApprox( hfTestObj->C0, 1e-5 ))
		printErr( (char*) "Inital MO matrix is not correct...", test, hfTestObj->C0 );

	cout << "Test 5 passed." << endl << endl;
}

/* --- Create and destroy new test constructs. --- */
sto_3g_water::sto_3g_water(double tol_e, double tol_dens)
{
	hfTestObj = new HartreeFock(tol_e, tol_dens);
}

void sto_3g_water::destroy(void)
{
	delete hfTestObj;
}

/* --- Test setup. --- */
int main(int argc, char* argv[])
{
	sto_3g_water tester( 0.01, 0.01 );
	tester.test1();
	tester.test2();
	tester.test3();
	tester.test4();
	tester.test5();
	cout << "Tests passed." << endl;
	return 0;
}

