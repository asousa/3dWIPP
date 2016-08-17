#include <Eigen/Core>

#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>
#include <iostream>

#include <consts.h>
#include <wipp.h>

#include <complex>
#include <cmath>

using namespace std;
using namespace Eigen;

extern "C" void sm_to_mag_d_(int* itime, double* x_in, double* x_out);
extern "C" void cart_to_pol_d_(double* x_in, double* lat, double* lon, double* radial);



// A port of the original raytracer's damping code, to compare to Forrest's 3d damping.

void damping_ngo(rayF ray) {
    int itime_in[2];
    double xin[3];
    double xout[3];
    double lat, lon, r;

    double lat_init, lon_init, r_init;
    double geom_fact;

    double B0, B0mag, kmag;
    double w;

    double wps2;
    complex<double> R, L, P, S, D; // Stix parameters

    complex<double> wcs;
    // complex<double> whs;
    double whs;
    VectorXd k, khat, Bhat;
    double kperp, kpar;
    double theta;




    itime_in[0] = 2012045;          // yyyyddd
    itime_in[1] = 1*(60*60*1000);   // time of day in msec     


    // // cout << "hello from damping";
    // cout << "file length is: " << ray.time.size() << "\n"; 

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




    // Loop over elements, do the calculation
    for (int ii=0; ii < ray.time.size(); ii++) {

        xin[0] = ray.pos[ii][0];
        xin[1] = ray.pos[ii][1];
        xin[2] = ray.pos[ii][2];
        
        // Map to magnetic dipole coordinates
        sm_to_mag_d_(itime_in, xin, xout);
        cart_to_pol_d_(xout, &lat, &lon, &r);
        // lat*= R2D;
        // lon*= R2D;
        r/= R_E;

        // printf("mag coords: %0.3f, %0.3f, %0.3f\n", R2D*lat, R2D*lon, r);


        // Calculate geometric factor (dipole coordinates)
        geom_fact = r_init * cos(lat_init) / (r * cos(lat));
        // printf("geometric factor: %0.3f\n",geom_fact);

        VectorXd B0  = Map<VectorXd>(ray.B0[ii].data(),3,1);
        VectorXd Ns  = Map<VectorXd>(ray.Ns[ii].data(),3,1);
        VectorXd n   = Map<VectorXd>(ray.n[ii].data(), 3,1);
        // MatrixXd Ns  = Map<MatrixXd>(ray.Ns[ii].data(),ray.nspec,1);
        // MatrixXd nus = Map<MatrixXd>(ray.nus[ii].data(),ray.nspec,1);
        // MatrixXd qs  = Map<MatrixXd>(ray.qs.data(),ray.nspec, 1);
        // MatrixXd ms  = Map<MatrixXd>(ray.ms.data(),ray.nspec, 1);

        // MatrixXd wps2 = Ns.array()*qs.array()*qs.array()/(ms.array()*EPS0);
        // MatrixXcd wcs  = ray.w*wps2.cwiseQuotient(ray.w + 1i*nus);
        // B0mag = l2_norm(ray.B0[ii]);
        B0mag = B0.norm();

        // ---------- Evaluate Stix parameters: ------
        wps2 = 0;
        wcs = 0;
        R = 1.;
        L = 1.;
        P = 1.;

        for (int jj=0; jj < ray.Ns[ii].size(); jj++) {
            
            // Ns*(Qs^2)/(ms*Eps_0)
            wps2 = ray.Ns[ii][jj]*pow(ray.qs[jj],2) \
                   /(ray.ms[jj]*EPS0);
            // qB/m
            whs  = ray.qs[jj]*B0mag/ray.ms[jj];

            wcs  = whs * ray.w/(ray.w + 1i*ray.nus[ii][jj]);

            R-= wps2/(ray.w*(ray.w + wcs));
            L-= wps2/(ray.w*(ray.w - wcs));
            P-= wps2/(ray.w*ray.w);
        }
        S = (R + L)/2.;
        D = (R - L)/2.;
        // --------------------------------------------


        // cout << B0.norm() << "\n";

        // Get K-vector and parallel / perp components WRT B0:
        k = n*ray.w/C;
        kmag = k.norm();

        Bhat = B0.array()/B0.norm();
        kpar = k.dot(Bhat); //k.array()*Bhat.array();

        kperp = (k - kpar*Bhat).norm();

        // Theta is the angle between parallel and perpendicular K
        theta = atan2(kperp, kpar);

        if (ray.time[ii] < 0.1) {
        cout << "Bhat: " << Bhat << "\n";
        cout << "K: " << k << "\n";
        cout << "kperp: " << kperp << " Kpar: " << kpar  << "\n";
        // cout << "K: " << k << "\n";
        printf("t: %0.2f w: %0.2f, B0: %e: Theta: %f\n",ray.time[ii],ray.w, B0mag,theta);
        printf("%0.2f %0.2f %0.2f %0.2f %0.2f\n",real(R),real(L),real(P),real(S),real(D));
    }

    }



} // damping_ngo




