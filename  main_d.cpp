#include <iostream>
#include <cmath>   /* sine function*/
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib> //atoi
#include "time.h" //

using namespace std;

double exact_sol( double x ) {return 1.0 - (1 - exp(-10))*x - exp(-10*x);}

int main(int argc, char* argv[]) {

    // Failtest for args
    if (argc <= 2) {
        cout << "Bad usage: " << argv[0] <<
             "read also output file on same line." << endl;
        exit(1);
    }

    // Declaring variables

    // Writing to file
    ofstream ofile;
    ofile.open((string(argv[1]) + "-" + string(argv[2]) + ".txt").c_str());

    // Declaration of variables
    int n = atoi(argv[2]); // second command line arg
    double *v_ = new double[n+2];
    double *u = new double[n+2];
    int spacing = 15;
    double *f = new double[n+2];
    double *x = new double[n+2];
    double *b = new double[n+2];
    double *f_tilde = new double[n+2];
    double *b_tilde = new double[n+2];
    double *z = new double[n+2];
    double h = 1.0/(n+1);

    //initialize elements
    for (int i = 0; i < n+1; i++) {
        x[i] = i*h;
        b[i] = 2.0;
        f[i] = h*h*100*exp(-10*x[i]);

    }

    // Calculating z and b_tilde
    for (int i = 1; i < n+2; i++) {
        b_tilde[i] = (i+1.0)/(i);
        z[i] = (i-1.0)/(i); // This is to reduce the amount of flops in f_tilde and v_
    }

    // Initial conditions
    f_tilde[1] = f[1];
    b_tilde[1] = b[1];

    // Looping over the algorithms from the methods-section - Forward substitution
    for (int i = 2; i < n+1; i++) {
        f_tilde[i] = f[i] + (z[i]*f_tilde[i-1]);
    }

    // 10th(and 100th and 1000th) element between the endpoints of the 12 points for n = 10 and so on.

    v_[n] = f_tilde[n]/b_tilde[n];

    // Backward substitution
    for (int i = n-1; i >= 1; i--) {
        v_[i] = (f_tilde[i] + v_[i+1])/b_tilde[i];
    }

    ofile << "x: " << setw(spacing + 5) << "v_: " << setw(spacing + 5)
          << "exact_sol: " << setw(spacing + 5) << "log10(Epsilon)" <<  endl;
    for (int i = 1; i < n; i++) {
         double Epsilon = fabs((v_[i] - exact_sol(x[i]))/exact_sol(x[i]));
             ofile << setw(15) << setprecision(8) << x[i];
             ofile << setw(15) << setprecision(8) << v_[i];
             ofile << setw(15) << setprecision(8) << exact_sol(x[i]);
             ofile << setw(15) << setprecision(8) << log10(Epsilon) << endl;
    }

    ofile.close();

    return 0;

};
