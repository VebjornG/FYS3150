#include <iostream>
#include <cmath>   /* sine function*/
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib> //atoi
#include "time.h" //
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {

    // Failtest for args
    if (argc <= 2) {
        cout << "Bad usage: " << argv[0] <<
             "read also output file on same line." << endl;
        exit(1);
    }

    // Declaring variables
    int n = atoi(argv[2]); // second command line arg
    int spacing = 15;

    vec b = zeros<vec>(n);
    vec f = zeros<vec>(n);
    vec x(n);
    mat A = zeros<mat>(n, n);
    double h = 1.0/(n+1);

    A(0, 0) = 2.0; A(0, 1) = -1.0; x(0) = h; b(0) = f(x(0));
    x(n-1) = x(0)+ (n-1)*h; b(n-1) = f(x(n-1));


    for (int i = 1; i < n-1; i++) {
       x(i) = x(i-1) + h;
        b(i) = f(i);
        A(i, i-1) = -1.0;
        A(i, i) = 2.0;
        A(i, i+1) = -1.0;
    }

    A(n-1,n-1) = 2.0; A(n-2, n-1) = -1.0; A(n-1, n-2) = -1.0;

    //Measuring runtime of LU decomposition
    clock_t start, finish;

    start = clock();


    // Solution for Av = f
    vec z = solve(A,f);

    finish = clock();

    // Time used calculating in the special case
    double runtime_LU = ((finish - start)/double(CLOCKS_PER_SEC));
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << setprecision(25) << setw(15) << "Runtime of the LU decomposition = " << runtime_LU << endl;

    // Writing to file
    ofstream ofile;
    ofile.open((string(argv[1]) + "-" + string(argv[2]) + ".txt").c_str());
    ofile << setw(spacing + 5) << setprecision(8) << "Time used running the LU decomposition: "<< runtime_LU << endl;
    ofile << setprecision(8) << A << endl;

    ofile.close();

    return 0;
};
