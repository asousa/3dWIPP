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


using namespace std;
using namespace Eigen;

void psd_model::initialize(char* inp_filename) {
    MATFile *pmat;
    DIR* dp;
    
    const mwSize *dims;
    
    // mxArray *MLT;
    // mxArray *J_perp;
    // mxArray *L;
    // mxArray *t_b;

    // dp = opendir(inp_dir_name);
    // int AE;
    // int E;
    // char *file;

    // ostringstream ss;
    
    // // Clear any previous data
    // MLT_v.clear();
    // J_perp_v.clear();
    // L_v.clear();
    // t_b_v.clear();

        // Load CRRES data (cleaned up pls!)
        pmat = matOpen(inp_filename, "r");
        if (pmat == NULL) {
            cout << "Error opening file at " << inp_filename << "\n";
        } else {
            crres_data = matGetVariable(pmat, "crres");
            num_entries = mxGetNumberOfElements(crres_data);

            cout << "Loaded " << inp_filename << "\n";

            // dims = mxGetDimensions(crres_data);

            // cout << "ndims: " << mxGetNumberOfFields(crres_data) << "\n";
            // cout << dims[0] << " " << dims[1] << "\n";
        }

        // // Load all the CRRES files in the directory
        // if (dp != NULL) {
        //     while (ep = readdir (dp)) {
        //         // cout << ep->d_name <<"\n";
        //         if (sscanf(ep->d_name, "crres_e%d_ae%d.mat",&E, &AE)) {
        //             // cout << meh << " " << AE << "\n";
        //             // cout << ep->d_name << "\n";
        //             ss.str("");
        //             ss.clear();
        //             ss << inp_dir_name << "/" << ep->d_name;
        //             // cout << ss.str() << "\n";

        //             pmat = matOpen(ss.str().c_str(), "r");
        //             if (pmat == NULL) {
        //                 printf("Error opening file %s\n", ep->d_name);
        //             } else {
        //                 // Load variables from this file:
        //                 MLT    = matGetVariable(pmat, "MLT");
        //                 J_perp = matGetVariable(pmat, "J_perp");
        //                 L      = matGetVariable(pmat, "L");
        //                 t_b    = matGetVariable(pmat, "t_b");

        //                 replace_NaNs(MLT);
        //                 replace_NaNs(J_perp);
        //                 replace_NaNs(L);
        //                 replace_NaNs(t_b);

        //                 // WELL: Now we have the shit, what do we do with it
        //                 MLT_v.push_back(MLT);
        //                 J_perp_v.push_back(J_perp);
        //                 L_v.push_back(L);
        //                 t_b_v.push_back(t_b);

        //                 // Energy and AE values from filename
        //                 energy_v.push_back((E*1e-3));
        //                 AE_v.push_back(AE);
        //             }
        //         } 
        //     }

        // dir_name = inp_dir_name;
        
        // } else {
        //     cout << "failed to read directory\n";
        // }
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



vector<double> psd_model::CRRES_fit_params(double L, double MLT, double AE_level) {
    // Loop through each of the CRRES entries, and find
    // the closest values for E and J_perp.
    // 
    // This is a port of "get_fit_params.m" by dgolden.
    
    int row, col;
    mxArray *L_mx;
    mxArray *MLT_mx;
    mxArray *Jp_mx;
    mxArray *E_mx;
    mxArray *AE_mx;


    double *MLTd;
    double *Ld;
    double *Jpd;
    MatrixXd MLT_data, L_data, Jp_data;
    double *E_d;
    double *AE_d;

    const mwSize *dims;

    // VectorXd log_Jp(num_entries);
    vector<double> log_Jp;
    vector<double> log_E;

    vector<double> fit_params;
    vector<double> p;   // polyfit coefficients
    double J0;
    double m;

    double AE_target = round(AE_level);
    cout << "AE_target: " << AE_target << "\n";

    // AE is either 1, 2, or 3 (after rounding)
    if ( (AE_target < 1) || (AE_target > 3) ) {
        cout << "AE_target out of range!\n";
        return fit_params;
    }

    for (int kk = 0; kk < num_entries; kk++) {
        // Get mxArrays from the structure
        MLT_mx = mxGetField(crres_data, kk, "MLT");         // MLT <mxn>
        L_mx   = mxGetField(crres_data, kk, "L");           // L-shell <mxn>
        Jp_mx  = mxGetField(crres_data, kk, "J_int");       // Jperp, interpolated <mxn>
        E_mx   = mxGetField(crres_data, kk, "E");           // Energy (keV) <double>
        AE_mx  = mxGetField(crres_data, kk, "AE");          // AE level     <double>

        // Get pointers to the sweet data within
        MLTd = mxGetPr(MLT_mx);
        Ld   = mxGetPr(L_mx);
        Jpd  = mxGetPr(Jp_mx);
        E_d  = mxGetPr(E_mx);
        AE_d = mxGetPr(AE_mx);

        if (AE_d[0] == AE_target) {

            // cout << "E: " << E_d[0] << " AE: "<< AE_d[0] << "\n";
            // get dimensions:
            dims = mxGetDimensions(MLT_mx);

            // Map to Eigen matrices
            MLT_data = Map<MatrixXd>(MLTd,dims[0], dims[1]);
            L_data   = Map<MatrixXd>(Ld,  dims[0], dims[1]);
            Jp_data  = Map<MatrixXd>(Jpd, dims[0], dims[1]);

            int number_of_dims=mxGetNumberOfDimensions(MLT_mx);
            int elements=mxGetNumberOfElements(MLT_mx);

            // Find the closest value to MLT in the first row of MLT_data:
            // (row, col indexes are probably constant for each run thru the CRRES data,
            // but let's leave this here for now unless it's obnoxiously slow)
            (MLT_data.row(0).array() - MLT).abs().minCoeff(&row);    

            // Same search for L vector:
            (L_data.col(0).array() - L).abs().minCoeff(&col);
           

            log_Jp.push_back(log(Jp_data(col, row)));
            log_E.push_back(log(E_d[0]));

        }   // AE level
    }   // Loop over each entry in crres_data
    
    cout << "Log_E: \n";
    print_vector(log_E);
    cout << "Log_Jperp: \n";
    print_vector(log_Jp);


    // Linear fit:
    polyfit(log_E, log_Jp, p, 1);
    
    cout << "Polyfit coeffs: \n";
    print_vector(p);

    J0 = exp(p[0]);
    m  = -1.0*p[1];

    cout << "J0: " << J0 << " m: " << m << "\n";
    n = 2*m + 2;
    // Watch out here! M_EL is our electron mass in consts.h,
    // but there's an "M_E" defined somewhere too.
    An = 2*J0/pow((0.5*(6.25e11)*M_EL), m-1);
    
    fit_params.push_back(n);
    fit_params.push_back(An);

    cout << "Fit Params: \n";
    print_vector(fit_params);

    return fit_params;
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


