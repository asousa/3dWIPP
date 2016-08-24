// wipp.h
#ifndef psd_H
#define psd_H

#include <vector>
#include <map>

#include <mat.h>

using namespace std;

// Plasmasphere phase-space density model (for damping)
class psd_model {
    char* dir_name;

    double An, n;

    // Vectors of pointers to 2d matrices for MLT, Jperp, L, and t_b
    // (from the CRRES data files)
    vector <mxArray*> MLT_v;
    vector <mxArray*> J_perp_v;
    vector <mxArray*> L_v;
    vector <mxArray*> t_b_v;
    vector <double> energy_v;
    vector <double> AE_v;

public: 
    pair<double,double> CRRES_fit_params(double L, double MLT, double AE_level);
    void initialize(char* input_dir);
    double suprathermal(double vperp, double vpar);
    double crres_psd(double vperp, double vpar);

private:
    void replace_NaNs(mxArray* arr);
    void count_NaNs(mxArray* arr);

};

#endif