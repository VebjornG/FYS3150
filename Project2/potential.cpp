#include "potential.h"
#include <cmath>

Potential::Potential(){
}

double Potential::HarmonicOscillator(double rho){
    return rho*rho;
}

double Potential::Coulomb(double w_r, double rho) {
    return w_r*w_r*rho*rho + 1./rho;
}

double Potential::HO_2e(double w_r, double rho) {
    return w_r*w_r*rho*rho;
}

