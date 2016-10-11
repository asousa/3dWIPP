// wipp.h
#ifndef wipp_H
#define wipp_H
#include <Eigen/Core>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>
#include <map>

#include <consts.h>

// #include <integrand.h>
// #include <psd_model.h>
// #include <mat.h>

using namespace std;

// Structure for holding an entire rayfile (yuge)
typedef struct rayF {
    int size;
    double w;                          // frequency (angular)
    double stopcond;                   // Stop condition
    double nspec;                      // Number of species in plasmasphere model
    // double ray_num;                 // Ray index number
    vector <double> time;            // Group time

    vector <vector <double> > pos;
    vector <vector <double> > vprel;
    vector <vector <double> > vgrel;
    vector <vector <double> > n;
    vector <vector <double> > B0;

    // Simulation time
    int iyr;
    int idoy;
    int isec;

    // Origin coordinates
    double in_radius, in_lat, in_lon, in_w;

    // Variable-length stuff (depending on number of constituents in model)
    vector <double> qs;    // species charge
    vector <double> ms;    // species mass
    vector <vector <double> > Ns;    // number density of species (m^-3)
    vector <vector <double> > nus;   // collision frequencies
    vector <double> damping;         // Damping vector (normalized to 1)

    double inp_pwr;                  // input power of ray

    // Stix parameters
    vector <double> stixP;
    vector <double> stixR;
    vector <double> stixL;
    vector <double> stixS;
    vector <double> stixD;
    vector <double> stixA;
    vector <double> stixB;    

} rayF;

// Structure for holding a single timestep of a ray (smol)
typedef struct rayT {
    double w;                          // frequency (angular)
    double nspec;                      // Number of species in plasmasphere model
    // double ray_num;                 // Ray index number
    double time;            // Group time

    double pos[3];
    double vprel[3];
    double vgrel[3];
    double n[3];
    double B0[3];

    // Variable-length stuff (depending on number of constituents in model)
    vector <double> qs;    // species charge
    vector <double> ms;    // species mass
    vector <double> Ns;    // number density of species (m^-3)
    vector <double> nus;   // collision frequencies
    
    double damping;        // Damping vector (normalized to 1)
    double inp_pwr;        // input power of ray

    // Stix parameters
    double stixP;
    double stixR;
    double stixL;
    double stixS;
    double stixD;
    double stixA;
    double stixB;

} rayT;


// rayfile loader
// map<int, rayF> read_rayfile(string fileName);

// Landau damping
// void damping_ngo(rayF &rayfile);

// Math functions
double l2_norm(vector<double> u);
vector<double> scalar_multiply(vector<double> u, double v);
double dot_product(vector<double>u, vector<double>v);
vector<double> add(vector<double>u, vector<double> v);

// Helpers
void print_vector(vector<double> u);
map<int, rayF> read_rayfile(string fileName);

// Science!
double input_power_scaling(double* flash_loc, double* ray_loc, double mag_lat, double w, double i0);
double ionoAbsorp(float lat, long f);
float interpPt(float *xI, float *yI, int n, float xO);

void interp_ray_fine(rayF** raylist, double n_x, double n_y, double n_z, int t_ind, rayT* out);
void calc_stix_parameters(rayF* ray);



// ---- Coordinate transforms ----
// lib_onera_desp (the Cospar IRBEM library, which is the internal support for SpacePy):

// extern "C" void coord_trans_vec1_(long* ntime, long* sysaxesIN, long* sysaxesOUT, long* iyr,long* idoy, double* secs, double* xIN, double* xOUT);
// // extern "C" void coord_trans1__(long* sysaxesIN, long* sysaxesOUT, long* iyr,long* idoy, long* secs, double* xIN, double* xOUT);

// extern "C" void sm2geo1_(long* iyr,long* idoy, double* secs, double* xSM, double* xGEO);
// extern "C" void geo2sm1_(long* iyr,long* idoy, double* secs, double* xGEO, double* xSM);

// extern "C" void geo2mag1_(long* iyr, double* xGEO, double* xMAG);
// // cartesian - spherical (trig terms in degrees!)
// extern "C" void car_sph_(double* xCAR, double* r, double* lat, double* loni);
// extern "C" void sph_car_(double* r, double* lat, double* loni, double* xCAR);


// ----- libxformd (Coordinate transform library used in the raytracer) -----
extern "C" void sm_to_geo_d_(int* itime, double* x_in, double* x_out);
extern "C" void geo_to_sm_d_(int* itime, double* x_in, double* x_out);

extern "C" void sm_to_mag_d_(int* itime, double* x_in, double* x_out);
extern "C" void mag_to_sm_d_(int* itime, double* x_in, double* x_out);

extern "C" void geo_to_mag_d_(int* itime, double* x_in, double* x_out);
extern "C" void mag_to_geo_d_(int* itime, double* x_in, double* x_out);

extern "C" void cart_to_pol_d_(double* x_in, double* lat, double* lon, double* radius);
extern "C" void pol_to_cart_d_(double* lat, double* lon, double* radius, double* x_out);





//  // Coordinate transforms:
// int sm2geo_init(char fname[],int *deg_ext,int *deg_day,int maxcoeff,double coeff[]);

// void sm2geo(double fday,int deg_want,                        
//              int deg_ext,int deg_day,double coeff[], 
//              double smcoeff[],double geocoeff[]);
// Porting the Damping Code:
// void damping_foust(rayF &rayfile);
// double integrand_wrapper(double x, void* data);

// double kp_to_pp(double kp);
// void polyfit(const vector<double> &xv, const vector<double> &yv, vector<double> &coeff, int order);

// // Plasmasphere phase-space density model (for damping)
// class psd_model {
//     char* dir_name;

//     double An, n;

//     // Vectors of pointers to 2d matrices for MLT, Jperp, L, and t_b
//     // (from the CRRES data files)
//     vector <mxArray*> MLT_v;
//     vector <mxArray*> J_perp_v;
//     vector <mxArray*> L_v;
//     vector <mxArray*> t_b_v;
//     vector <double> energy_v;
//     vector <double> AE_v;

// public: 
//     pair<double,double> CRRES_fit_params(double L, double MLT, double AE_level);
//     void initialize(char* input_dir);
//     double suprathermal(double vperp, double vpar);
//     double crres_psd(double vperp, double vpar);

// private:
//     void replace_NaNs(mxArray* arr);
//     void count_NaNs(mxArray* arr);

// };


#endif
