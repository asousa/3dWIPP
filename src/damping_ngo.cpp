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

#include <complex>
#include <cmath>

using namespace std;
using namespace Eigen;

extern "C" void sm_to_mag_d_(int* itime, double* x_in, double* x_out);
extern "C" void cart_to_pol_d_(double* x_in, double* lat, double* lon, double* radial);



// Resonance mode (I think)
#define     M_RES 0
// I think these are fit params for the suprathermal electron distribution
#define     P_DIST      0.0
#define     Q_DIST      2.0
#define     AN_CM_DIST  2E5
#define     V0_DIST     1.0


// A port of the original raytracer's damping code, to compare to Forrest's 3d damping.

void damping_ngo(rayF &ray) {
    int itime_in[2];
    double xin[3];
    double xout[3];
    double lat, lon, r;

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


    // I don't understand what this is yet.
    AN = AN_CM_DIST * pow( 10.0 , (12.0-(4.0*Q_DIST)) );


    // Initialize ray power with zeros
    // VectorXd ray_pwr = VectorXd::Zero(ray.time.size());
    VectorXd ray_pwr = VectorXd::Zero(ray.time.size());

    // Initial power
    // ray.damping.push_back(1.0);
    ray.damping = vector<double> (ray.time.size(), 0.0);
    ray.damping[0] = 1.0;

    // Loop over elements, do the calculation
    for (int ii=1; ii < ray.time.size(); ii++) {

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

        B0    = Map<VectorXd>(ray.B0[ii].data(),3,1);
        Ns    = Map<VectorXd>(ray.Ns[ii].data(),3,1);
        n_vec = Map<VectorXd>(ray.n[ii].data(), 3,1);

        B0mag = B0.norm();

        // ---------- Evaluate Stix parameters: ------
        wps2 = 0;
        // wcs = 0;
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
        // --------------------------------------------


        // cout << B0.norm() << "\n";

        // Get K-vector and parallel / perp components WRT B0:
        k = n_vec*ray.w/C;
        Bhat = B0.array()/B0.norm();
        kpar = k.dot(Bhat); //k.array()*Bhat.array();
        kperp = (k - kpar*Bhat).norm();

        // Theta is the angle between parallel and perpendicular K
        theta = atan2(kperp, kpar);

        // Some trig.
        sin_th = sin(theta);
        cos_th = cos(theta);
        sin_th_sq = pow(sin_th,2);
        cos_th_sq = pow(cos_th,2);

        a = S*sin_th_sq + P*cos_th_sq;
        b = R*L*sin_th_sq + P*S*(1+cos_th_sq);

        dens = Ns[0];  // electron number density (first species is electrons, right?)
                     

        DELTA = 1./dens;
        Omega = (Q_EL/M_EL)*B0mag;  // Electron gyrofrequency, radians (qb/m)
            
        // refractive index 
        // (equation and norm of the vector match! Sticking with the vector)
        // n = sqrt(  pow(C/ray.w,2)*(pow(kperp,2) + pow(kpar,2))  );
        n = n_vec.norm();
        n_sq = pow(n, 2);

        // Plasma frequency
        wp = Q_EL*sqrt( dens / (EPS0*M_EL) );  
        
        const1 = DELTA/ (4*n*(2*a*n_sq - b));
        Vm = (ray.w + M_RES*Omega) / kpar;  
        Vm_sq = pow(Vm,2);

        c1 = (4* AN * pow(PI,2) * pow(wp,2) ) / ( ray.w * kpar ) ;

        chi1 = 0.0;

        // Here's where we integrate over v_perp!
        for(v_perp=10; v_perp<C ; v_perp+=v_step ) {

            v_perp_sq = pow(v_perp,2);
            bessel_arg = kperp * v_perp / Omega ;
            bessel_term = (R-n_sq) * jn((M_RES-1),bessel_arg) + 
                (L-n_sq) * jn((M_RES+1),bessel_arg) ;
            v_tot_sq = v_perp_sq + Vm_sq + V0_DIST ; 
            distrib_term = P_DIST*v_tot_sq  -  v_perp_sq*Q_DIST; 
            velocity_term = pow( v_tot_sq , (Q_DIST+1) );


            //TERM1:
            c2 = (n_sq*sin_th_sq - P) / (2.0*(S - n_sq) );
            term1 =   c2 * c1 * 
                pow( bessel_term , 2 ) *
                pow( v_perp , twoPPlusOne ) / (ray.w*velocity_term) *
                ( ray.w*distrib_term - kpar*Vm*P_DIST*v_tot_sq ) ;

            //TERM2:
            g2 = (1+ M_RES * Omega / ray.w) * (Q_DIST * Vm * v_perp) + 
                 (M_RES*Omega*Vm / (ray.w * v_perp) ) * distrib_term;

            term2 = 2.*c1 * ( (S-n_sq*cos_th_sq)*(S-n_sq) - pow(D,2) ) *
                pow( jn(M_RES,bessel_arg) , 2 ) * Vm * pow(v_perp_sq,P_DIST)/
                velocity_term * g2 ;

            //TERM3:
            c3 = 2.*n_sq*sin_th*cos_th ;
            term3 = c3 * c1 * pow(v_perp,twoPPlusOne)*jn(M_RES,bessel_arg) /
            velocity_term * bessel_term * g2 ;

            // chi1 += term1 + term2 + term3;
            // Original WIPP code doesn't put the abs() here, but
            // this makes the damping match. Hmm... hmm.
            chi1 += abs(term1 + term2 + term3) ;
        }   // v_para


        // ds -- distance traveled by ray in this timestep (in meters?)
        pos      = Map<VectorXd>(ray.pos[ii].data(),3,1);
        pos_prev = Map<VectorXd>(ray.pos[ii-1].data(), 3,1);
        
        ds = (pos - pos_prev).norm();

        k_im = (chi1 * const1 * v_step * ray.w / C);
        damp += k_im * ds;
        // damp += (k_im * ds);

        // ray_pwr[ii] = geom_fact * exp(2.0 * damp) * init_pwr;
        // ray.damping.push_back(exp(2.0 * damp) * ray.damping[0]);
        ray.damping[ii] = (geom_fact * exp(2.0 * damp) * ray.damping[0]);


        // printf("t: %0.2f  pwr: %f\n",ray.time[ii], ray.damping[ii]);

        // if (ray.time[ii] < 0.1) {
        //     cout << "n: " << n << ", n_vec mag: " << n_vec.norm() << "\n";
        //     cout << "Bhat: " << Bhat << "\n";
        //     cout << "K: " << k << "\n";
        //     cout << "kperp: " << kperp << " Kpar: " << kpar  << "\n";
        //     cout << "ds: " << ds << "\n";
        //     // cout << "K: " << k << "\n";
        //     printf("t: %0.2f w: %0.2f, B0: %e: Theta: %f\n",ray.time[ii],ray.w, B0mag,theta);
        //     printf("%0.2f %0.2f %0.2f %0.2f %0.2f\n", R, L, P, S, D);
        // }



    }
    // ray.damping = ray_pwr.data();


} // damping_ngo




