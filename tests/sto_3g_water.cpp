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
    void test6(void);
    void test7(void);
    void test8(void);
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
	cout << "Running test5: Checking Initial Fock Matrix..." << endl;
	MatrixXd test(7, 7);
	readMatrix2d("data/fock_t.dat", &test);
	reflectMatrix2d(&test, 1);
	if(!test.isApprox( hfTestObj->F0, 1e-5 ))
		printErr( (char*) "Inital fock matrix is not correct...", test, hfTestObj->F0 );

	cout << "Running test5: Checking Initial MO Matrix..." << endl;
	readMatrix2d("data/mo_i.dat", &test);
	if(!test.isApprox( hfTestObj->C0, 1e-5 ))
		printErr( (char*) "Inital MO matrix is not correct...", test, hfTestObj->C0 );

	cout << "Running test5: Checking Initial Density Matrix..." << endl;
	readMatrix2d("data/dens_i.dat", &test);
	reflectMatrix2d(&test, 1);
	if(!test.isApprox( hfTestObj->D0, 1e-5 ))
		printErr( (char*) "Inital MO matrix is not correct...", test, hfTestObj->D0 );

	cout << "Test 5 passed." << endl << endl;
}

//Test 6 - Testing New Fock Matrix
void sto_3g_water::test6(void)
{
	cout << "Running Test 6." << endl;
	cout << "Running test6: Checking New Fock Matrix..." << endl;
	MatrixXd test(7, 7);
	readMatrix2d("data/fock_n.dat", &test);
	reflectMatrix2d(&test, 1);
	if(!test.isApprox( hfTestObj->F0, 1e-5 ))
		printErr( (char*) "New Fock matrix is not correct...", test, hfTestObj->F0 );
	cout << "Test 6 passed. " << endl << endl;
}

//Test 7 - Testing Final Energies
void sto_3g_water::test7(void)
{
	cout << "Running Test 7." << endl;
	cout << "Running test7: Checking Final energies..." << endl;
	double e_elec = -82.944446990003;
	double e_tot = -74.942079928192;
	if( abs( hfTestObj->eelec - e_elec ) > 1e-5 )
	{
		cout << "Initial electron energy is incorrect, should be " << e_elec << ", but yours is " << hfTestObj->eelec << "." << endl;
		assert(false);
	}

	if( abs( hfTestObj->etot - e_tot ) > 1e-5 )
	{
		cout << "Initial total energy is incorrect, should be " << e_tot << ", but yours is " << hfTestObj->etot << "." << endl;
		assert(false);
	}
	cout << "Test 7 passed. " << endl << endl;
}

//Test 8 - Testing Final Values
void sto_3g_water::test8(void)
{
	cout << "Running Test 8." << endl;
	cout << "Running test8: Checking final values..." << endl;

	double mu_x = 0.000000000000;
	double mu_y = 0.603521296525;
	double mu_z = mu_x;
	double dm_tot = 0.603521296525;
	double q0 = -0.253146052405;
	double q1 = 0.126573026202;
	double q2 = q1;

	if( abs( hfTestObj->mu_x - mu_x ) > 1e-5 )
	{
		cout << "Mu_x is incorrect, should be " << mu_x << ", but yours is " << hfTestObj->mu_x << "." << endl;
		assert(false);
	}

	if( abs( hfTestObj->mu_y - mu_y ) > 1e-5 )
	{
		cout << "Mu_y is incorrect, should be " << mu_y << ", but yours is " << hfTestObj->mu_y << "." << endl;
		assert(false);
	}

	if( abs( hfTestObj->mu_z - mu_z ) > 1e-5 )
	{
		cout << "Mu_z is incorrect, should be " << mu_z << ", but yours is " << hfTestObj->mu_z << "." << endl;
		assert(false);
	}

	if( abs( hfTestObj->tot_dip_moment - dm_tot ) > 1e-5 )
	{
		cout << "Total dipole moment is incorrect, should be " << dm_tot << ", but yours is " << hfTestObj->tot_dip_moment << "." << endl;
		assert(false);
	}

	if( abs( hfTestObj->q0 - q0 ) > 1e-5 )
	{
		cout << "Atom 0 charge is incorrect, should be " << q0 << ", but yours is " << hfTestObj->q0 << "." << endl;
		assert(false);
	}

	if( abs( hfTestObj->q1 - q1 ) > 1e-5 )
	{
		cout << "Atom 1 charge is incorrect, should be " << q1 << ", but yours is " << hfTestObj->q1 << "." << endl;
		assert(false);
	}

	if( abs( hfTestObj->q2 - q2 ) > 1e-5 )
	{
		cout << "Atom 2 charge is incorrect, should be " << q2 << ", but yours is " << hfTestObj->q2 << "." << endl;
		assert(false);
	}

	cout << "Test 8 passed. " << endl << endl;
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
	tester.test6();
	tester.test7();
	tester.test8();

	cout << "Tests passed." << endl;
	return 0;
}
