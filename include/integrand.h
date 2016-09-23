// integrand.h
#ifndef integrand_H
#define integrand_H

#include <vector>
#include <map>

#include <mat.h>


using namespace std;

// Integrand (the core integration in damping_foust.cpp)
class integrand {

double kperp, kpar, w, wch, qs, Ns, ms, nus, B0;
int m_low, m_hi;
double theta, n, n_sq, ct, st, st_sq;
double R, L, P, S;

psd_model f;

public: 
    void initialize(psd_model& f_in, double kperp_in, double kpar_in, double w_in, 
                   int m_low_in, int m_hi_in, double wch_in,
                   double R_in, double L_in, double P_in, double S_in);

    // simple evaluation
    double evaluate(double vperp);
    // normalized by speed of light
    double evaluate_vperpnorm(double vperp);
    // remapped onto finite range [0, 1]
    double evaluate_t(double vperp, void* data);

private:
    double fG1(double vperp, double vpar);
    double fG2(double vperp, double vpar, double m);


};

#endif
