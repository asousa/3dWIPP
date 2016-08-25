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

#include <complex>
#include <cmath>

using namespace std;
using namespace Eigen;

extern "C" void sm_to_mag_d_(int* itime, double* x_in, double* x_out);
extern "C" void cart_to_pol_d_(double* x_in, double* lat, double* lon, double* radial);

class psd_model; 


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

    double lat_init, lon_init, r_init;
    double geom_fact;

    double B0mag, kmag;
    // double w;

    double wps2;
    double R, L, P, S, D; // Stix parameters
    double a, b;

    // complex<double> wcs;
    // complex<double> whs;
    double whs;
    VectorXd k, khat, Bhat;
    double kperp, kpar;
    double theta;

    double sin_th, cos_th, sin_th_sq, cos_th_sq;
    double dens;
    double DELTA, Omega, n, n_sq, wp;
    double Vm, Vm_sq, c1;
    double AN;
    double chi1;
    double v_perp, v_step, n_steps;

    double v_perp_sq, bessel_arg, v_tot_sq;
    double bessel_term;

    double distrib_term, velocity_term;
    double twoPPlusOne =2*P_DIST + 1;

    double const1;

    double c2;
    double term1;
    double g2, term2;
    double c3, term3;

    double ds;
    Vector3d pos;
    Vector3d pos_prev;
    Vector3d B0;
    VectorXd Ns;
    VectorXd n_vec;

    double k_im, damp;
    double init_pwr = 1.0;   // Initial ray power


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
    //                     fe = @(vperp, vpar) crres_polar_hybrid_psd(vperp, vpar, n_fit, An_fit, L, L_pp);    char crres_data_file[100] = "damping_standalone/crres_data_clean/crres_clean.mat";
    
    psd_model psd;
    psd.initialize(crres_data_file);
    //               [n,  An]                        L,  MLT, AE
    vector<double> fit_params = psd.CRRES_fit_params(2.2, 1, 3);
    // Inputs are: 
    //                     vperp, vpar, n,             An,            L, L_pp
    double f = psd.hybrid_psd(1,  2,    fit_params[0], fit_params[1], 3, 2);
    cout << "f: " << f << "\n";
    // -----------------------------------------------------------------------


    // Here's where you left off!
    //  --- damp_single, line 136:  wce_h and onward.


} // damping_ngo




