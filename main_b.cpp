#include <iostream>
#include <cmath>   /* sine function*/
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib> //atoi
#include "time.h" //

using namespace std;


// General solution with 8n flops

int main(int argc, char* argv[]) {

    // Failtest for args
    if (argc <= 2) {
        cout << "Bad usage: " << argv[0] <<
             "read also output file on same line." << endl;
        exit(1);
    }

    // Declaring variables
    int n = atoi(argv[2]); // second command line arg
    double *v = new double[n+2];
    int spacing = 15;
    double *a = new double[n+2];
    double *f = new double[n+2];
    double *x = new double[n+2];
    double *b = new double[n+2];
    double *c = new double[n+2];
    double *f_tilde = new double[n+2];
    double *b_tilde = new double[n+2];
    //double *runtime_gen = new double[n];
    double h = 1.0/(n+1);

    //Setting up the matrices
    for (int i = 0; i < n+2; i++) {
        x[i] = i*h;
        v[i] = 0.0;
        a[i] = -1.0;
        c[i] = -1.0;
        b[i] = 2.0;
        f[i] = h*h*100*exp(-10*x[i]);
    }
    //Initial conditions
    f_tilde[1] = f[1];
    b_tilde[1] = b[1];

    // Forward substitution
    for (int i = 2; i < n+1; i++) {
        b_tilde[i] = b[i] - (a[i-1]*c[i-1]/b_tilde[i-1]);
        f_tilde[i] = f[i] - (a[i-1]*f_tilde[i-1]/b_tilde[i-1]);
    }

    // 10th(and 100th and 1000th) element between the endpoints of the 12 points for n = 10 and so on.
    v[n] = f_tilde[n]/b_tilde[n];

      // Backward substitution
    for (int i = n-1; i >= 1; i--) {
        v[i] = (f_tilde[i] - v[i+1]*c[i])/b_tilde[i];
    }

    // Writing to file
    ofstream ofile;
    ofile.open((string(argv[1]) + "-" + string(argv[2]) + ".txt").c_str());
    ofile << "          x:" << setw(spacing + 5) << "          v:" << "          f_tilde:" << endl;
    for (int i = 0; i < n+2; i++) {
        ofile << setw(spacing + 5) << setprecision(8) << v[i];
        ofile << setw(spacing + 5) << setprecision(8) << x[i];
        ofile << setw(spacing + 5) << setprecision(8) <<  f_tilde[i] << endl;
    }

    ofile.close();

    return 0;
};
