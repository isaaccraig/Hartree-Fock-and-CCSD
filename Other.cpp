//The mass header file

#include "mass.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <armadillo>
#include <complex>


/**************************************************************/
//Building the initial guess density matrix
/**************************************************************/

	arma::mat I_Gs=S_sqrt.t()*H_c*S_sqrt;
	I_Gs.print("\nThis is the initial Fock matrix(orthogonal basis)\n");
//Note that this contains the inital orbital range
	arma::vec eig11;
	arma::mat eigenvec;
	arma::eig_sym(eig11,eigenvec,I_Gs);
	real(eig11).print("\nThe eigenvalues (initially)\n");
	real(eigenvec).print("\nThe eigenvectors (initially)\n");
	//eigvec1.print("\nThe initial MO coefficient Matrix is \n");
//Transforming the eigenvectors into the original basis
	arma::mat eig_AO = S_sqrt*eigenvec;
	real(eig_AO).print("\nThe initial MO cofficients\n");
//The density matrix
//Since only the first MO of water are occupied so the lowest five
//MO are taken into consideration and they are effectively squared.
	eig_AO.shed_cols(5,6);
	arma::mat den_I = eig_AO*eig_AO.t();
	real(den_I).print("\nThe initial guess density matrix\n");

/**************************************************************/
//Computing the new Fock Matrix and starting the SCF loop
/**************************************************************/
	arma::mat density_initial=den_I;
	double energy_last=energy;
//setting the precision here
	cout.precision(15);
	arma::mat F_mv=H_c;
        arma::mat coeff_car=den_I;
//the SCF loop starts here
	for(int ii=1; ii<100000; ii++)
	{
		F_mv=H_c;
		for(int i=0;i<H_c.n_rows;i++)
		{
			for(int j=0; j<H_c.n_cols; j++)
			{
				for(int k=0;k<H_c.n_cols;k++)
				{
					for(int l=0;l<H_c.n_cols;l++)
					{
						F_mv(i,j)+=den_I(k,l)*(2*d_I(index_call(i+1,j+1,k+1,l+1))-d_I(index_call(i+1,k+1,j+1,l+1)));
					}
				}
			}
		}
		//F_mv.print("\nThe Fock matrix with 2 electron integrals \n");
		arma::mat F_store=S_sqrt.t()*F_mv*S_sqrt;
		arma::vec eig_store;
		arma::mat eigenvec_store;
		arma::eig_sym(eig_store,eigenvec_store,F_store);
		arma:: mat eig_AO_store=S_sqrt*eigenvec_store;
		coeff_car=eig_AO_store;
		//eig_AO_store.print("\n");
		eig_AO_store.shed_cols(5,6);
		arma:: mat density_store=eig_AO_store*eig_AO_store.t();
		double energy1=0.0;
		for(int i=0; i< H_c.n_rows; i++)
		{
			for(int j=0; j<H_c.n_cols; j++)
			{
				energy1+=density_store(i,j)*(H_c(i,j)+F_mv(i,j));
			}
		}
		double d_ene=energy_last-energy1;
		double rms=0.0;
		double rms_last=0.0;
		for(int i=0; i<density_store.n_rows;i++)
		{
			for (int j=0;j<density_store.n_cols;j++)
			{
				rms+=pow(density_store(i,j)-density_initial(i,j),2);
			}
		}
		rms=pow(rms,0.5);
		den_I=density_store;
		printf("\n %d\t  %10.12f \t %10.12f \t %10.12f \t %10.12f \n ",ii,energy1,energy1+enuc,d_ene,rms);
		//F_mv.print("\n");
		//den_I.print("\nThe density matrix after quitting the loop\n");
		if (abs(d_ene)<1e-12 && abs(rms_last-rms)< 1e-12 )
			break;
//Changing the stuff here
		density_initial=density_store;
		energy_last=energy1;
		rms_last=rms;

	}
/**************************************************************/
//The MO-basis Fock-Fatrix
/**************************************************************/
	//F_mv.print("\n");
	//coeff_car.print("\n");
	arma::mat F_MO2=coeff_car.t()*F_mv*coeff_car;
	cout<<"\n\nChecking for the Diagonal form of Fork Matrix in the basis of the Molecular Orbitals\n\n";

	for (int i=0; i<F_MO2.n_rows;i++)
	{
		for(int j=0; j<F_MO2.n_cols;j++)
		{
			printf("%15.7f",F_MO2(i,j));
		}
		cout<<endl;
	}

/********************************************************************/
//One Electron Properties- electronic contribution to dipole operator
/*********************************************************************/
//Reading the files
	arma::mat mux_file;
	mux_file.load("mux.dat");
	arma::mat muy_file;
	muy_file.load("muy.dat");
	arma::mat muz_file;
	muz_file.load("muz.dat");
//Constructing the matrix (one dimensional based on index)
	int ele=muz_file(muz_file.n_rows-1,0);
	arma::mat mux=arma::zeros(ele,ele);
	arma::mat muy=arma::zeros(ele,ele);
	arma::mat muz=arma::zeros(ele,ele);
	for(int i=0; i< mux_file.n_rows;i++)
	{
		int x1=mux_file(i,0);
		int x2=mux_file(i,1);
		double val_x=mux_file(i,2);
		mux(x1-1,x2-1)=val_x;

		int y1=muy_file(i,0);
		int y2=muy_file(i,1);
		double val_y=muy_file(i,2);
		muy(y1-1,y2-1)=val_y;

		int z1=muz_file(i,0);
		int z2=muz_file(i,1);
		double val_z=muz_file(i,2);
		muz(z1-1,z2-1)=val_z;

	}
	mux=symmatl(mux);
	muy=symmatl(muy);
	muz=symmatl(muz);
	mux.print("\nThis is the mu_x integrals\n");
	muy.print("\nThis is the mu_y integrals\n");
	muz.print("\nThis is the mu_z integrals\n");

	double mu_x=0.0;
	double mu_y=0.0;
	double mu_z=0.0;
	double N_x=0.0;
	double N_y=0.0;
	double N_z=0.0;

	for(int i=0; i<mux.n_rows; i++)
	{
		for(int j=0; j<mux.n_cols;j++)
		{
			mu_x+=2*den_I(i,j)*mux(i,j);
			mu_y+=2*den_I(i,j)*muy(i,j);
			mu_z+=2*den_I(i,j)*muz(i,j);
		}
	}
	for(int i=0; i<n_atoms; i++)
	{
		N_x+=z_val(i)*geom(i,0);
		N_y+=z_val(i)*geom(i,1);
		N_z+=z_val(i)*geom(i,2);
	}
	cout.precision(15);
	printf(" \n %s \t %8.5f \n %s \t %8.5f \n %s \t %8.5f \n","The dipole moment operators are: \n\n mu_x",mu_x+N_x,"mu_y",mu_y+N_y,"mu_z",mu_z+N_z);
	printf("\n %s \t %8.5f\n", "The total dipole moment is :", mu_x+N_x+mu_y+N_y+mu_z+N_z);

}



void geometry(arma::mat &geom, int &n_atoms, arma::vec &z_val)
{
	ifstream myreadfile1("geom.dat");
	if(myreadfile1.is_open())
		{
			cout<<"\nReading the geometry file now\n"<<endl;
			myreadfile1 >> n_atoms;
			geom=arma::zeros(n_atoms,n_atoms);
			z_val=arma::zeros(n_atoms);
			for(int i=0; i<n_atoms; i++)
			{
				myreadfile1 >> z_val(i) >> geom(i,0) >> geom(i,1) >>
				geom(i,2);
			}
		}
}


int index1(int i, int j)
{
	int ij=i*(i+1)/2+j;
	return ij;
}


int index_call(int i, int j , int k, int l)
{
	int ij,kl,ijkl;
	if(i>j)
	{
		 ij=index1(i,j);
	}
	else
	{
		 ij=index1(j,i);
	}
	if (k>l)
	{
		 kl=index1(k,l);
	}
	else
	{
		kl=index1(l,k);
	}
	if(ij > kl)
	{
		 ijkl=ij*(ij+1)/2+(kl);
	}
	else //if (ij < kl)

	{
		 ijkl=kl*(kl+1)/2+(ij);
	}

	return ijkl;
}


© 2018 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
API
Training
Shop
Blog
About
