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

#include <sys/types.h>
#include <dirent.h>
#include <regex.h>

#include "mat.h"



void psd_model::initialize(char* inp_dir_name) {
    MATFile *pmat;
    DIR* dp;
    struct dirent *ep;     
    const char *name;
    mxArray *MLT;
    mxArray *J_perp;
    mxArray *L;
    mxArray *t_b;

    dp = opendir(inp_dir_name);
    int AE;
    int E;
    char *file;

    ostringstream ss;
    
    // Clear any previous data
    MLT_v.clear();
    J_perp_v.clear();
    L_v.clear();
    t_b_v.clear();


        // Load all the CRRES files in the directory
        if (dp != NULL) {
            while (ep = readdir (dp)) {
                // cout << ep->d_name <<"\n";
                if (sscanf(ep->d_name, "crres_e%d_ae%d.mat",&E, &AE)) {
                    // cout << meh << " " << AE << "\n";
                    // cout << ep->d_name << "\n";
                    ss.str("");
                    ss.clear();
                    ss << inp_dir_name << "/" << ep->d_name;
                    // cout << ss.str() << "\n";

                    pmat = matOpen(ss.str().c_str(), "r");
                    if (pmat == NULL) {
                        printf("Error opening file %s\n", ep->d_name);
                    } else {
                        // Load variables from this file:
                        MLT    = matGetVariable(pmat, "MLT");
                        J_perp = matGetVariable(pmat, "J_perp");
                        L      = matGetVariable(pmat, "L");
                        t_b    = matGetVariable(pmat, "t_b");

                        replace_NaNs(MLT);
                        replace_NaNs(J_perp);
                        replace_NaNs(L);
                        replace_NaNs(t_b);

                        // WELL: Now we have the shit, what do we do with it
                        MLT_v.push_back(MLT);
                        J_perp_v.push_back(J_perp);
                        L_v.push_back(L);
                        t_b_v.push_back(t_b);

                        // Energy and AE values from filename
                        energy_v.push_back((E*1e-3));
                        AE_v.push_back(AE);
                    }
                } 
            }

        dir_name = inp_dir_name;
        
        } else {
            cout << "failed to read directory\n";
        }
}


void psd_model::replace_NaNs(mxArray* arr) {
    // Allegedly the CRRES data uses 1e12 to denote a NaN.
    /* Declare variables */ 
    size_t elements;
    mwSize j,cmplx;
    mwSize number_of_dims;
    mwSize nnz=0, count=0; 
    double *pr, *pi, *pind;
    const mwSize *dim_array;

    /* Get the data */
    pr=(double *)mxGetPr(arr);
    pi=(double *)mxGetPi(arr);
    cmplx = ((pi==NULL) ? 0 : 1);

    number_of_dims=mxGetNumberOfDimensions(arr);
    elements=mxGetNumberOfElements(arr);
    // cout << " n_dims: " << number_of_dims;
    // cout << " elements: " << elements << "\n";

    for(int j=0;j<elements;j++) {
        if (pr[j]==1e12) {
            pr[j] = NAN;
            // cout << " nan @ " << j;
        }
    }
}

void psd_model::count_NaNs(mxArray* arr) {
    // Allegedly the CRRES data uses 1e12 to denote a NaN.
    /* Declare variables */ 
    size_t elements;
    mwSize j,cmplx;
    mwSize number_of_dims;
    mwSize nnz=0, count=0; 
    double *pr, *pi, *pind;
    const mwSize *dim_array;

    /* Get the data */
    pr=(double *)mxGetPr(arr);
    pi=(double *)mxGetPi(arr);
    cmplx = ((pi==NULL) ? 0 : 1);

    // cout << "complex? " << cmplx;

    number_of_dims=mxGetNumberOfDimensions(arr);
    elements=mxGetNumberOfElements(arr);
    // cout << " n_dims: " << number_of_dims;
    // cout << " elements: " << elements << "\n";
    int nanTotal = 0;
    for(int j=0;j<elements;j++) {
        if (mxIsNaN(pr[j])) {
            nanTotal++;
        }
    }
    cout << "Total NaNs: " << nanTotal << "\n";
}



pair<double,double> psd_model::CRRES_fit_params(double L, double MLT, double AE_level) {

    cout << "number of files: " << MLT_v.size() << "\n";

    print_vector(energy_v);
    print_vector(AE_v);
}



double psd_model::suprathermal(double vperp, double vpar) {
    // Suprathermal distribution.
    // This is what we'll use inside the plasmasphere.

    double v, f;

    const double a = 4.9e5;
    const double b = 8.3e14;
    const double c = 5.4e23;
    // % Just a crutch to avoid the singularity
    const double v0=1;

    // % Convert to cm/s
    v = 100.*sqrt(vperp*vperp + vpar*vpar + v0);
    f = (a/pow(v,4) -b/pow(v,5) + c/pow(v,6));

    // % Convert to s^3/m^6 from s^3/cm^6 
    f = f*pow(100.,6);
    return f;
}

double psd_model::crres_psd(double vperp, double vpar) {
    // Suprathermal distribution, model-driven.
    // This is what we'll use outside the plasmasphere.

    // % Get phase space density using CRRES suprathermal fluxes from Bortnik
    // % [2007]
    // % 
    // % vperp and vpar in m/s

    // % By Daniel Golden (dgolden1 at stanford dot edu) May 2010
    // % $Id: crres_psd.m 1205 2011-03-04 23:15:40Z dgolden $

    double v, f;

    // % Just a crutch to avoid the singularity
    const double v0=1;

    // % Convert to cm/s
    v = 100.*sqrt(vperp*vperp + vpar*vpar + v0);

    f = An/pow(v, n);

    // % Convert to s^3/m^6 from s^3/cm^6 
    f = f*pow(100.,6);
    return f;
}


