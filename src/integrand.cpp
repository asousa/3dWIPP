#include <Eigen/Core>

#include <stdio.h>
#include <math.h>

#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>


#include <consts.h>
#include <wipp.h>
#include <psd_model.h>
#include <integrand.h>

#include <limits.h>

using namespace std;
using namespace Eigen;

void integrand::initialize(psd_model& f_in, double kperp_in, double kpar_in, 
                           double w_in, double n_in, int m_low_in, int m_hi_in, double wch_in,
                           double R_in, double L_in, double P_in, double S_in) {

    // Current parameters:
    kperp = kperp_in;
    kpar = kpar_in;
    w = w_in;
    m_low = m_low_in;   // Lowest resonance mode  (0 is landau, +- 1 is cyclotron)
    m_hi  = m_hi_in;    // Highest resonance mode
    wch = wch_in;   

    // Pointer to psd model:
    f = f_in;

    // Stix params:
    R = R_in;
    L = L_in;
    P = P_in;
    S = S_in;  

// % angle of wavenormal with respect to B0
    theta = atan2(kperp,kpar);

    // % Refractive index
    // n_sq = ( C*C/(w*w) ) * (kperp*kperp + kpar*kpar);
    // n = sqrt(n_sq);
    n = n_in;
    n_sq = n*n;
    // cout << "n (in): " << n << "\n";

    // % cos(theta)
    ct = cos(theta);
    // % sin(theta)
    st = sin(theta);

    st_sq = st*st;
    // cout << "theta: " << theta << "\n";
}


double integrand::evaluate(double vperp) {
// evaluate the integral at vperp:

double suminteg = 0;
double vpar;
double Jm, Jm_m1, Jm_p1;
double G1, G2;
double ret;

    
    // Evaluate for each resonance between m_low and m_hi:
    // Jm = jn(1, 2);
    // cout << "Jm: " << Jm << "\n";

    for (int m = m_low; m < m_hi; m++) {
        // evaluate Bessel functions
        Jm = jn(m, kperp*vperp/wch);
        Jm_m1 = jn(m + 1, kperp*vperp/wch);
        Jm_p1 = jn(m - 1, kperp*vperp/wch);

        vpar = (w - m*wch)/kpar;

        G1 = fG1(vperp, vpar);
        G2 = fG2(vperp, vpar, m);

    // % Chen's derivation.  Kennel has accidentally swapped (R-n^2) and (L-n^2)!
    suminteg += (G1* ( (P-n_sq*st_sq)*(2*(L-n_sq)*vperp*Jm_p1*Jm_p1 + 2*vperp*(R-n_sq)*Jm_m1*Jm_m1 +
                                       n_sq*st_sq*vperp*pow(Jm_p1-Jm_m1, 2)) 
                       -n_sq*ct*st*(2*vpar*Jm*(Jm_p1*(R-n_sq)+Jm_m1*(L-n_sq)) +
                                    n_sq*ct*st*vperp*pow(Jm_p1-Jm_m1, 2)))
         +G2*(4*vpar*Jm*((L-n_sq)*(R-n_sq)+n_sq*st_sq*(S-n_sq))
              -2*n_sq*ct*st*((R-n_sq)*vperp*Jm_m1+(L-n_sq)*vperp*Jm_p1)));
    }

    ret = -2*PI*PI*((Q_EL*Q_EL/M_EL/EPS0)/(w*fabs(kpar))) * suminteg*vperp;
    return ret;
}


double integrand::evaluate_vperpnorm(double vperp) {
    return C*evaluate(vperp*C);
}

double integrand::evaluate_t(double t) {
    double eps = DBL_EPSILON; // double-precision epsilon, from limits.h
    // cout << "evaluate_t: " << t << "\n";
    return ((1+eps)/(t*t + eps))*evaluate_vperpnorm((1-t+eps)/(t+eps));
}



double integrand::fG1(double vperp, double vpar) {
    const double DEL=1e-8;
    double eps = DBL_EPSILON;  // double-precision epsilon, from limits.h
    double G1;
    // cout << "eps: " << eps << "\n";
    double d = DEL*abs(vperp);
    
    if( d < 10*eps ) {
      d = 10*eps;
    }

    double df_dvperp = ( f.hybrid_psd(vperp+d, vpar) -
                         f.hybrid_psd(vperp-d, vpar) ) 
                        / (2*d);

    d = DEL*abs(vpar);

    if ( d < 10*eps ) {
        d = 10*eps;
    }

    double df_dvpar = ( f.hybrid_psd(vperp, vpar+d) -
                        f.hybrid_psd(vperp, vpar-d) ) 
                        / (2*d);
    
    G1 = df_dvperp - (kpar/w)*(vpar*df_dvperp - vperp*df_dvpar);
    // cout << "df_perp: " << df_dvperp << " df_par: " << df_dvpar << " G1: " << G1 << "\n";
    return G1;
}

double integrand::fG2(double vperp, double vpar, double m) {
    const double DEL=1e-8;
    double eps = DBL_EPSILON;  // double-precision epsilon, from limits.h
    double G1;
    // cout << "eps: " << eps << "\n";
    double d = DEL*abs(vperp);
    
    if( d < 10*eps ) {
      d = 10*eps;
    }

    double df_dvperp = ( f.hybrid_psd(vperp+d, vpar) -
                         f.hybrid_psd(vperp-d, vpar) ) 
                        / (2*d);

    d = DEL*abs(vpar);

    if ( d < 10*eps ) {
        d = 10*eps;
    }

    double df_dvpar = ( f.hybrid_psd(vperp, vpar+d) -
                        f.hybrid_psd(vperp, vpar-d) ) 
                        / (2*d);
    
    double Jm = jn(m, kperp*vperp/wch);

    double G2 = Jm*(df_dvpar-(m*wch+eps)/(w*vperp+eps)*(vpar*df_dvperp - vperp*df_dvpar));

    return G2;

}


