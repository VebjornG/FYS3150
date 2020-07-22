#include <iostream>
#include <cmath>   /* sine function*/
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib> //atoi
#include "time.h" //

using namespace std;


// General solution with 8n flops
double gen_solv(int n, double* v) {

    // Declaration
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

    clock_t start, finish;

    // start clock for gen_solv
    start = clock();

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

    // end clock
    finish = clock();

    // Time used calculating in the general case
    /*double runtime_gen = 0;
    for (int i = 0; i <= n; i++) {
        double runtime_gen[i] = ((finish - start)/double(CLOCKS_PER_SEC));
    }*/

    //double average_gen = runtime_gen/n;
    double runtime_gen = ((finish - start)/double(CLOCKS_PER_SEC));
    //cout << average_gen << endl;
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << setprecision(25) << setw(15) << "The time used in the general case = " << runtime_gen << endl;
    return runtime_gen;
};
// End gen_sol

// Start special solver
double Special_solv(int n, double* v_) {

    // Declaration
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


    clock_t start, finish;

    // start clock for Special_solv
    start = clock();

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

    // end clock
    finish = clock();

    // Time used calculating in the special case
    double runtime_special = ((finish - start)/double(CLOCKS_PER_SEC));
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << setprecision(25) << setw(15) << "The time used in the special case = " << runtime_special << endl;
    return runtime_special;
};


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
    double *v_ = new double[n+2];
    int spacing = 15;

    // Calling the solvers to write the times and v, v_ to file
    double t_general = gen_solv(n, v);
    double t_special = Special_solv(n, v_);

    // Writing to file
    ofstream ofile;
    ofile.open((string(argv[1]) + "-" + string(argv[2]) + ".txt").c_str());
    ofile << setw(spacing + 5) << setprecision(8) << "Time used gen_solv: " << t_general << endl;
    ofile << setw(spacing + 5) << setprecision(8) << "Time used Special_solv: "<< t_special << endl;
    ofile << "v:" << setw(spacing + 5) << "v_:" << endl;
    for (int i = 0; i < n+2; i++) {
        ofile << setw(spacing) << setprecision(8) << v[i];
        ofile << setw(spacing + 5) << setprecision(8) << v_[i] << endl;
    }

    ofile.close();

    return 0;
};