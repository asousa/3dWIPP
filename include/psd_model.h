// wipp.h
#ifndef psd_H
#define psd_H

#include <vector>
#include <map>

#include <mat.h>

using namespace std;

// Plasmasphere phase-space density model (for damping)
class psd_model {
    char* filename;

    double An, n;

    double L_sh, L_pp;  // L-shell and L-plasmapause
    double MLT, AE_level;
    
    mxArray* crres_data;
    int num_entries;

    // vector<double> fit_params;
    double n_fit, An_fit;   // Fit params

    // Vectors of pointers to 2d matrices for MLT, Jperp, L, and t_b
    // (from the CRRES data files)
    // vector <mxArray*> MLT_v;
    // vector <mxArray*> J_perp_v;
    // vector <mxArray*> L_v;
    // vector <mxArray*> t_b_v;
    // vector <double> energy_v;
    // vector <double> AE_v;

public: 
    vector<double> CRRES_fit_params(double L, double MLT, double AE_level);
    void initialize(char* inp_filename);
    double suprathermal(double vperp, double vpar);
    double crres_psd(double vperp, double vpar, double n, double An);
    // double hybrid_psd(double vperp, double vpar, double n_fit, double An_fit, double L, double L_pp);
    double hybrid_psd(double vperp, double vpar);

    void set_params(double L_sh_in, double L_pp_in, double MLT_in, double AE_level_in);


private:
    void replace_NaNs(mxArray* arr);
    void count_NaNs(mxArray* arr);

};

#endif
