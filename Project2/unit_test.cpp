//#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#ifdef CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <armadillo>
#include "systemsolver.h"
#include <iostream>

using namespace arma;
using namespace std;


TEST_CASE( "Find max value", "Approximation" ) {
    int N = 3;

    // test matrix
    mat test_mat;
    test_mat     << 1 << 3 << 2 << endr
                 << 3 << 4 << 6 << endr
                 << 2 << 6 << 7 << endr;

    // getting max value from maxOffDiag
    SystemSolver *system = new SystemSolver(N);
    double max_value = system->Offdiag_maxvalue(test_mat);

    int &k = system->k;
    int &l = system->l;

    // Checking whether the Offdiag_maxvalue function
    // cout << test_mat(1,2) << endl;
    REQUIRE(max_value == Approx(6));
    cout << "k = " << k << ", " << "l = " << l << endl;
    cout << test_mat << endl;

    // Hardcoding the elements to zero and testing for the next values
    test_mat(k,l) = 0; test_mat(l,k) = 0;
    max_value = system->Offdiag_maxvalue(test_mat);
    REQUIRE(max_value == Approx(3));
    cout << test_mat << endl;

    test_mat(k,l) = 0; test_mat(l,k) = 0;
    max_value = system->Offdiag_maxvalue(test_mat);
    REQUIRE(max_value == Approx(2));
    cout << test_mat << endl;

    test_mat(k,l) = 0; test_mat(l,k) = 0;
    max_value = system->Offdiag_maxvalue(test_mat);
    REQUIRE(max_value == Approx(0));
    cout << test_mat << endl;
}

TEST_CASE( "Testing eigenvalues", "Approximation" ) {
    int N = 3;

    // test matrix
    mat test_mat;
    test_mat << 1 << 5 << 6 << endr
             << 5 << 5 << 1 << endr
             << 6 << 1 << 2 << endr;

    mat test_mat_Z = eye<mat>(N,N);

    double lambda_0 = -5.50093442;
    double lambda_1 = 2.71301099;
    double lambda_2 = 10.78792343;

    SystemSolver* test_matrix = new SystemSolver(N);

    // Naming the matrices made in the class SystemSolver

    test_matrix->jacobi_method(test_mat, test_mat_Z); // Diagonal matrix
    vec eigvals = zeros<vec>(N);

    for (int i = 0; i < N; i++) {
        eigvals(i) = test_mat(i, i);
    }
    eigvals = sort(eigvals);
    for (int i = 0; i < 3; i++) {
        cout << "Lamba_" << i << "= " << eigvals(i) << endl;
    }
    REQUIRE(eigvals(0) == Approx(lambda_0));
    REQUIRE(eigvals(1) == Approx(lambda_1));
    REQUIRE(eigvals(2) == Approx(lambda_2));
}
#endif
