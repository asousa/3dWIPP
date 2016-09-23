#include <Eigen/Core>

#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>


#include <consts.h>
#include <wipp.h>
#include <psd_model.h>
#include <integrand.h>
#include <gauss_legendre.h>

#include <complex>
#include <cmath>

using namespace std;
using namespace Eigen;

extern "C" void sm_to_mag_d_(int* itime, double* x_in, double* x_out);
extern "C" void cart_to_pol_d_(double* x_in, double* lat, double* lon, double* radial);

class psd_model; 
class integrand;

// Resonance mode (I think)
#define     M_RES 0
// I think these are fit params for the suprathermal electron distribution
#define     P_DIST      0.0
#define     Q_DIST      2.0
#define     AN_CM_DIST  2E5
#define     V0_DIST     1.0

#define     KP          4.0  // Add this as an input or something, you doof.

// A port of Forrest's 3d damping code.
void damping_foust(rayF &ray) {
    int itime_in[2];
    double xin[3];
    double xout[3];
    double lat, lon, r;
    double L_sh;
    double MLT;

    double lat_init, lon_init, r_init;
    double AN;
    double v_step, n_steps;

    Vector3d pos;
    Vector3d B0;
    Vector3d n_vec;
    Vector3d k;
    Vector3d Ns;

    Vector3d Bhat;

    double kpar, kperp;
    double theta;
    double sin_th, cos_th, sin_th_sq, cos_th_sq;

    double n;

    double k_im, damp;
    double init_pwr = 1.0;   // Initial ray power

    double fs;
    double wce_h;
    double kmag;

    double R, L, P, S, D; // Stix parameters
    double a, b;
    double B0mag;
    double wps2, whs;
    // Change this to an input, you doof
    double AE_level = 3;

    itime_in[0] = 2012045;          // yyyyddd
    itime_in[1] = 1*(60*60*1000);   // time of day in msec     

    n_steps = 500.0;
    v_step = C/n_steps; //<- change back to this!

    // Get ray launch location in mag dipole coordinates:
    xin[0] = ray.pos[0][0];
    xin[1] = ray.pos[0][1];
    xin[2] = ray.pos[0][2];
    
    // // Map to magnetic dipole coordinates
    sm_to_mag_d_(itime_in, xin, xout);
    cart_to_pol_d_(xout, &lat_init, &lon_init, &r_init);
    // printf("pos size: %d\n",ray.pos[0].size());
    r_init/= R_E;

    printf("mag coords: %0.3f, %0.3f, %0.3f\n", R2D*lat_init, R2D*lon_init, r_init);


    // I don't understand what this is yet.
    AN = AN_CM_DIST * pow( 10.0 , (12.0-(4.0*Q_DIST)) );


    // Initialize ray power with zeros
    VectorXd ray_pwr = VectorXd::Zero(ray.time.size());

    // Initial power
    ray.damping = vector<double> (ray.time.size(), 0.0);
    ray.damping[0] = 1.0;

    // get L-shell of plasmapause, based on kp:
    double L_pp = kp_to_pp(KP);

    cout << L_pp << "\n";


    // ----------- Here's how we get the phase-space density model: ----------
    //  (equivalent to:   [n_fit, An_fit] = get_fit_params(L, MLT, AE_level, false);
    //                     fe = @(vperp, vpar) crres_polar_hybrid_psd(vperp, vpar, n_fit, An_fit, L, L_pp);    
    char crres_data_file[100] = "damping_standalone/crres_data_clean/crres_clean.mat";
    
    psd_model psd;
    psd.initialize(crres_data_file);

    // Step through the ray:
    for (int ii=0; ii < ray.time.size(); ii++) {

        B0    = Map<VectorXd>(ray.B0[ii].data(), 3,1);
        pos   = Map<VectorXd>(ray.pos[ii].data(),3,1);
        n_vec = Map<VectorXd>(ray.n[ii].data(),3,1);
        Ns    = Map<VectorXd>(ray.Ns[ii].data(),3,1);

        // Get local L-shell:
        lat = atan(pos[2]/sqrt(pow(pos[0],2) + pow(pos[1],2)));
        r = pos.norm();
        L_sh = r/pow(cos(lat),2)/R_E;

        // Get MLT:
        MLT = fmod((atan2(pos[1], pos[0]) + PI)/(2*PI)*24, 24); //MLT in hours; 0 MLT is in -x direction

        // Set the current location parameters for the density model:
        psd.set_params(L_sh, L_pp, MLT, AE_level);

        wce_h = Q_EL*B0.norm()/M_EL;

        k = n_vec*ray.w/C;
        kmag = k.norm();
        Bhat = B0.array()/B0.norm();
        kpar = k.dot(Bhat); //k.array()*Bhat.array();
        kperp = (k - kpar*Bhat).norm();

        B0mag = B0.norm();

        // Print some shit:
        cout << "t: " << ray.time[ii] << " MLT: " << MLT << " lat: " << R2D*lat << " L_sh: " << L_sh << "\n";
        
        // ------- spatialdamping.m -----------
        // Theta is the angle between parallel and perpendicular K
        theta = atan2(kperp, kpar);

        // Some trig.
        sin_th = sin(theta);
        cos_th = cos(theta);
        sin_th_sq = pow(sin_th,2);
        cos_th_sq = pow(cos_th,2);


        n = n_vec.norm();

        // ---------- Evaluate Stix parameters: ------
        wps2 = 0;
        R = 1.;
        L = 1.;
        P = 1.;

        for (int jj=0; jj < ray.Ns[ii].size(); jj++) {
            
            // Ns*(Qs^2)/(ms*Eps_0)
            wps2 = ray.Ns[ii][jj]*pow(ray.qs[jj],2) \
                   /(ray.ms[jj]*EPS0);
            // qB/m
            whs  = ray.qs[jj]*B0mag/ray.ms[jj];

            // Complex modification to whs -- for now, ignore. (8.19.2016)
            // wcs  = whs * ray.w/(ray.w + 1i*ray.nus[ii][jj]);

            R-= wps2/(ray.w*(ray.w + whs));
            L-= wps2/(ray.w*(ray.w - whs));
            P-= wps2/(ray.w*ray.w);
        }
        S = (R + L)/2.;
        D = (R - L)/2.;

        a = S*sin_th_sq + P*cos_th_sq;
        b = R*L*sin_th_sq + P*S*(1+cos_th_sq);
        // --------------------------------------------

        // ---------- hot_dispersion_imag.m ------

        integrand integ;
        int m_low = -1;
        int m_hi  = 1;

        integ.initialize(psd, kperp, kpar, 
                        ray.w, m_low, m_hi, wce_h, 
                        R, L, P, S);

        // cout << integ.evaluate_t(0.5) << "\n";

        // // Integrate it!
        // double tmp_integ = gauss_legendre(2, eff, NULL, 0, 2*PI);
        // cout << "tmp integ: " << tmp_integ << "\n";


    } // Step through ray


for (double t=0; t < 2*PI; t+=0.01) {
    double tmp_integ = gauss_legendre(20, eff, NULL, 0, t);
    cout << "t: " << t << " integ: " << tmp_integ << "\n";

}

} // damping_ngo


// Here's a working place to write an integrand. But how do we point it to
// a method of an object? Hmm.
double eff(double x, void* data) {
    return x;
}
