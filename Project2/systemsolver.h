#ifndef SYSTEMSOLVER_H
#define SYSTEMSOLVER_H
#include <armadillo>
#include <string>

class SystemSolver
{
private:
    int N;
public:
    int k, l;
    SystemSolver(int size);
    arma::mat init(double rho_max, std::string potential);
    void jacobi_method(arma::mat &A, arma::mat &Z);
    double Offdiag_maxvalue(arma::mat A);
    void rotation(arma::mat &A, arma::mat &Z);
    void eig_symmetric(arma::mat &A);
};
#endif // SYSTEMSOLVER_H
