#ifndef CATCH_CONFIG_MAIN
#include "potential.h"
#include "systemsolver.h"
#include <time.h>
#include <iostream>
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iomanip>

using namespace arma;
using namespace std;

int main(){
    double rho_max = 30;
    int n = 350;
    SystemSolver* system = new SystemSolver(n);
    string potential = "CO";


    // Naming the matrices made in the class SystemSolver
    mat A = system->init(rho_max, potential);
    mat Z = eye<mat>(n,n); // Diagonal matrix
    system->jacobi_method(A, Z);

    // Armadillo funciton eig_sym
    system->eig_symmetric(A);

    vec eigvals = zeros<vec>(n);

    // Filling in the vector eigvals with the eigenvalues
    for (int i = 0; i < n; i++) {
        eigvals(i) = A(i, i);
    }

    // Finding the right indices for the eigenvalues corresponding
    // to the right eigenvectors
    uvec indices = sort_index(eigvals);

    // Sorting the eigenvalues
    eigvals = sort(eigvals);
    for (int i = 0; i < 3; i++) {
        cout << "Lamba_" << i << "= " << eigvals(i) << endl;
    }


    // Writiting the eigenvalues to file
    ofstream ofile;
    ofile.open("Eigenvals.txt");
    for (int i = 0; i < n; i++) {
        ofile << setprecision(8) << eigvals(i) << endl;
    }
    ofile.close();

    // Writiting the three first eigenvectors to file
    ofile.open("Eigenvectors_" + potential + "_" + "wr=0.25" + ".txt");
    for (int i = 0; i < n; i++) {
            ofile << setw(20) << setprecision(8) << Z(i, indices(0));
            ofile << setw(20) << setprecision(8) << Z(i, indices(1));
            ofile << setw(20) << setprecision(8) << Z(i, indices(2)) << endl;
    }
    ofile.close();

    return 0;
}
#endif
