#include "systemsolver.h"
#include "potential.h"
#include <armadillo>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <time.h>

using namespace arma;
using namespace std;

SystemSolver::SystemSolver(int size) {
    N = size;
}

// Constructing the matrix A
mat SystemSolver::init(double rho_max, string potential) {

    Potential* pot = new Potential();
    mat A = zeros<mat>(N, N);
    double h = rho_max/N;
    double e_i = 1.0/(h*h);
    for (int i = 0; i < N; i++){
        // Checking whether we're working with the interacting or non-interacting case
        if (potential == "HO") {
            A(i, i) = 2.*e_i + pot->HarmonicOscillator((i+1)*h);
        }
        else if (potential == "CO") {
            A(i, i) = 2.*e_i + pot->Coulomb(0.25, (i+1)*h);
        }
        else if (potential == "HO_2e") {
            A(i, i) = 2.*e_i + pot->HO_2e(0.25, (i+1)*h);
        }
        if (i < N-1) {
            A(i+1, i) = -e_i;
            A(i, i+1) = -e_i;
        }
    }
    return A;
}

void SystemSolver::jacobi_method(mat &A, mat &Z) {

    double epsilon = 1.0e-8;
    double max_num_iterations = (double) N * (double) N * (double) N;
    double Offdiag_max = Offdiag_maxvalue(A);
    int iterations = 0;

    // Runtime for the Jacobi method
    clock_t start, finish;
    start = clock();

    // Making the offdiagonal elements smaller than epsilon
    while (fabs(Offdiag_max) > epsilon && (double) iterations < max_num_iterations) {
        Offdiag_max = Offdiag_maxvalue(A);
        rotation(A, Z);
        iterations++;
    }
    finish = clock();
    std::cout << "Number of iterations: " << iterations << "\n";

    // Time used
    double runtime_JM = ((finish - start)/double(CLOCKS_PER_SEC));
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << setprecision(9) << setw(15) << "Runtime of Jacobi method = " << runtime_JM << endl;
}


// Function to find the maximum offdiagonal value
double SystemSolver::Offdiag_maxvalue(mat A) {
    double maxvalue = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            if (fabs(A(i, j)) > maxvalue) {
                maxvalue = fabs(A(i, j));
                this->k = j;
                this->l = i;
            }
        }
    }
    return maxvalue;
}

// Rotation matrix to find the values of cos and sin
void SystemSolver::rotation(mat &A, mat &Z) {
    double s, c;
    if (A(k, l) != 0.0) {
        double t, tau;
        tau = (A(l, l) - A(k, k))/(2.*A(k, l));
        if (tau > 0) {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }
        c = 1.0/sqrt(1.0 + t*t);
        s = c*t;
    }
    else {
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k, k);
    a_ll = A(l, l);

    //Changing the matrix elements with indices k and l
    A(k, k) = c*c*a_kk - 2.0*c*s*A(k, l) + s*s*a_ll;
    A(l, l) = s*s*a_kk + 2.0*c*s*A(k, l) + c*c*a_ll;
    A(k, l) = 0.0;
    A(l, k) = 0.0;

    //Changing the remaining elements
    for (int i = 0; i < N; i++) {
        if (i != k && i != l) {
            a_ik = A(i, k);
            a_il = A(i, l);
            A(i, k) = c*a_ik - s*a_il;
            A(k, i) = A(i, k);
            A(i, l) = c*a_il + s*a_ik;
            A(l, i) = A(i, l);
        }

        // Computing the new eigenvectors
        r_ik = Z(i, k);
        r_il = Z(i, l);
        Z(i, k) = c*r_ik - s*r_il;
        Z(i, l) = c*r_il + s*r_ik;
    }
}
void SystemSolver::eig_symmetric(mat &A) {

    // Calculate computation time and eigenvectors with standard Armadillo library:
    clock_t start_arma, finish_arma;
    start_arma = clock();
    mat eigvec;
    vec eigval;
    eig_sym(eigval, eigvec, A);
    finish_arma = clock();

    // Time used
    double runtime_arma = (finish_arma - start_arma)/(double)CLOCKS_PER_SEC;
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << setprecision(15) << setw(15) << "Runtime of eig_sym = " << runtime_arma << endl;
}
