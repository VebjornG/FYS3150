/*
    Problem 4b).

    Program to compute the mean energy, magnetization,
    specific heat capacity and susceptibility as functions of T
    on a 2x2 lattice.
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <random>
#include "lib.h"

using namespace std;
ofstream ofile;

// Inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}

// Function to read in data from screen
void read_input(int&, double&);

// Function to initialise energy and magnetization
void initialize(int, double, int **, double&, double&);

// The metropolis algorithm
void Metropolis(int, long&, int **, double&, double&, double *);

// prints to file the results of the calculations
void output(int, int, double, double *);

int main(int argc, char* argv[]) {
    char *outfilename;
    long idum;
    int **spin_matrix, n_spins, mcs;
    double w[17], average[5], temperature, E, M;

    // Read in output file, abort if there are too few command-line arguments
    if( argc <= 1 ){
      cout << "Bad Usage: " << argv[0] <<
        " read also output file on same line" << endl;
      exit(1);
    }
    else{
      outfilename=argv[1];
    }
    ofile.open(outfilename);

    //    Read in initial values such as size of lattice, temp and cycles
    read_input(n_spins, temperature);
    idum = -1; // random starting point
    int n[] = {1E4, 1E5, 1E6, 1E7};
    for (int j=0; j <= 3; j++) {
      spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
      mcs = n[j];

      //    initialise energy and magnetization
      E = M = 0.;

      // setup array for possible energy changes
      for( int de =-8; de <= 8; de++) w[de+8] = 0;
      for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);

      // initialise array for expectation values
      for( int i = 0; i < 5; i++) average[i] = 0.;
      initialize(n_spins, temperature, spin_matrix, E, M);

      // start Monte Carlo computation
      for (int cycles = 1; cycles <= mcs; cycles++){
      Metropolis(n_spins, idum, spin_matrix, E, M, w);

      // Update the expectation values
      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M;  average[4] += fabs(M);
      }

    // Print results
    output(n_spins, mcs, temperature, average);
    free_matrix((void **) spin_matrix);
    }
    ofile.close();
    return 0;
}

// read in input data
void read_input(int& n_spins, double& temperature)
{
  cout << "Lattice size (number of spins where x and y equal): ";
  cin >> n_spins;
  cout << "Temperature with dimension energy: ";
  cin >> temperature;
}

// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, double temperature, int **spin_matrix,
        double& E, double& M){

  // setup spin matrix and intial magnetization
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      spin_matrix[y][x] = 1; // spin orientation for the ground state
      M +=  (double) spin_matrix[y][x];
    }
  }
  // setup initial energy
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      E -=  (double) spin_matrix[y][x]*
    (spin_matrix[periodic(y,n_spins,-1)][x] +
     spin_matrix[y][periodic(x,n_spins,-1)]);
    }
  }
}

void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w){

  // loop over all spins
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      int ix = (int) (ran1(&idum)*(double)n_spins);
      int iy = (int) (ran1(&idum)*(double)n_spins);
      int deltaE =  2*spin_matrix[iy][ix]*
                    (spin_matrix[iy][periodic(ix,n_spins,-1)]+
                    spin_matrix[periodic(iy,n_spins,-1)][ix] +
                    spin_matrix[iy][periodic(ix,n_spins,1)] +
                    spin_matrix[periodic(iy,n_spins,1)][ix]);
      if ( ran1(&idum) <= w[deltaE+8] ) {
        spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
        M += (double) 2*spin_matrix[iy][ix];
        E += (double) deltaE;
      }
    }
  }
} // end of Metropolis sampling over spins

void output(int n_spins, int mcs, double temperature, double *average){

  double norm = 1/((double) (mcs));  // divided by total number of cycles
  double norm2 = 1.0/(n_spins*n_spins);   // divide by the total number of spins
  double T = temperature;
  double Eaverage = average[0]*norm;
  double E2average = average[1]*norm;
  double Maverage = average[2]*norm;
  double M2average = average[3]*norm;
  double Mabsaverage = average[4]*norm;

  // Exact results for T = 1.0 (J = 1):
  double cosh_fac = (3.0+cosh(8/T));
  double E_exact = -8.0*sinh(8/T)/cosh_fac*norm2;
  double absM_exact = (2.0*exp(8.0/T)+4.0)/cosh_fac*norm2;
  double Cv_exact = (8.0/T)*(8.0/T)*(1.0+3*cosh(8/T))/(cosh_fac*cosh_fac)*norm2;
  double Suscept_exact = 1/T*(12.0+8.0*exp(8.0/T)+8*cosh(8.0/T))/(cosh_fac*cosh_fac)*norm2;

  // all expectation values are per spin, divide by 1/n_spins/n_spins
  double Evariance = (E2average - Eaverage*Eaverage)/n_spins/n_spins;
  double Mvariance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "Number of Monte Carlo trials: " << mcs << endl;
  ofile << "Temperature:      Energy:      Heat capacity:  Magnetization:  Susceptibility:  abs(Magnetization):" << endl;
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << Eaverage/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
  ofile << setw(15) << setprecision(8) << Maverage/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Mvariance/temperature;
  ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins << endl;

  // Relative errors:
  ofile << "Diff(E):          Diff(C_V):    Diff(Chi):     Diff(abs(M)):" << endl;
  ofile << setw(15) << setprecision(8) << abs(Eaverage*norm2-E_exact)/abs(E_exact);
  ofile << setw(15) << setprecision(8) << abs(Evariance/temperature/temperature-Cv_exact)/Cv_exact;
  ofile << setw(15) << setprecision(8) << abs(Mvariance/temperature-Suscept_exact)/Suscept_exact;
  ofile << setw(15) << setprecision(8) << abs(Mabsaverage*norm2-absM_exact)/absM_exact << endl;
} // end output function
