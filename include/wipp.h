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
#include <bmodel.h>

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

// Structure for holding parameters for a single EA segment 
// (planes perpendicular to field line of interest, where we
// will calculate scattering at)
typedef struct EA_segment {

    Eigen::Vector3d ea_norm;           // Vector normal to crossing plane
    Eigen::Vector3d field_line_pos;    // Location in plane where field line intersects

    double L_sh;                       // L shell
    double radius;                     // Radius around field line to consider a crossing

    double dist_to_n;
    double dist_to_s;
    // Should we do Stix parameters here? They're a constant of the background
    // medium, not of the wave... but we'd have to get the plasmasphere model.

} EA_segment;


// rayfile loader
// map<int, rayF> read_rayfile(string fileName);

// Landau damping
// void damping_ngo(rayF &rayfile);

// Math functions
double l2_norm(vector<double> u);
double norm(double u[], int size);
vector<double> scalar_multiply(vector<double> u, double v);
double dot_product(vector<double>u, vector<double>v);
vector<double> add(vector<double>u, vector<double> v);

// Helpers
void print_vector(vector<double> u);
void print_array(double* arr, double len);

map<int, rayF> read_rayfile(string fileName);

// Science!
double input_power_scaling(double* flash_loc, double* ray_loc, double mag_lat, double w, double i0);
double ionoAbsorp(float lat, long f);
float interpPt(float *xI, float *yI, int n, float xO);

void interp_ray_fine(rayF** raylist, double n_x, double n_y, double n_z, int t_ind, rayT* out);
void calc_stix_parameters(rayF* ray);
void init_EA_array(EA_segment* EA_array, double lat, double lon, int iyr, int idoy, double isec);


// ---- Coordinate transforms ----
// lib_onera_desp (the Cospar IRBEM library, which is the internal support for SpacePy):

// extern "C" void coord_trans_vec1_(long* ntime, long* sysaxesIN, long* sysaxesOUT, long* iyr,long* idoy, double* secs, double* xIN, double* xOUT);
// // extern "C" void coord_trans1__(long* sysaxesIN, long* sysaxesOUT, long* iyr,long* idoy, long* secs, double* xIN, double* xOUT);

// extern "C" void sm2geo1_(long* iyr,long* idoy, double* secs, double* xSM, double* xGEO);
// extern "C" void geo2sm1_(long* iyr,long* idoy, double* secs, double* xGEO, double* xSM);

extern "C" void geo2mag1_(int* iyr, double* xGEO, double* xMAG);

// // cartesian - spherical (trig terms in degrees!)
extern "C" void car_sph_(double* xCAR, double* r, double* lat, double* loni);
extern "C" void sph_car_(double* r, double* lat, double* loni, double* xCAR);


// ----- libxformd (Coordinate transform library used in the raytracer) -----
extern "C" void sm_to_geo_d_(int* itime, double* x_in, double* x_out);
extern "C" void geo_to_sm_d_(int* itime, double* x_in, double* x_out);

extern "C" void sm_to_mag_d_(int* itime, double* x_in, double* x_out);
extern "C" void mag_to_sm_d_(int* itime, double* x_in, double* x_out);

extern "C" void geo_to_mag_d_(int* itime, double* x_in, double* x_out);
extern "C" void mag_to_geo_d_(int* itime, double* x_in, double* x_out);

// (in radians)
extern "C" void cart_to_pol_d_(double* x_in, double* lat, double* lon, double* radius);
extern "C" void pol_to_cart_d_(double* lat, double* lon, double* radius, double* x_out);

extern "C" void sm_to_gsm_d_(int* itime, double* x_in, double* x_out);
extern "C" void gsm_to_sm_d_(int* itime, double* x_in, double* x_out);

// ---- My own transforms ----
// In-place cartesian / polar transforms. 
void carsph(double x[3]); 
void sphcar(double x[3]); 
void cardeg(double x[3]);
void degcar(double x[3]);

// In-place mapping of a data field between cartesian / polar frames.
void transform_data_sph2car(double lat, double lon, double d_in[3], double d_out[3]);
void transform_data_car2sph(double lat, double lon, double d_in[3], double d_out[3]);
void transform_data_geo2mag(int itime_in[2], double d_in[3], double d_out[3]);
void transform_data_mag2geo(int itime_in[2], double d_in[3], double d_out[3]);

// ----- liboneradesp (Cospar IRBEM library -- IGRF wrapper, field line tracer, etc)
extern "C" void get_field_multi_(int* ntime, int* kext, int options[5],
                                 int* sysaxes, int* iyear, int* idoy, double*ut,
                                 double* x1, double* x2, double* x3,
                                 double maginput[25], double Bgeo[], double* Bl);



// // ----- libgeopackd (Tyganenko's transform library + IGRF12)
extern "C" void igrf_geo_08_(double* r, double* theta, double* phi, double* Br, double* Btheta, double* Bphi);
extern "C" void igrf_gsw_08_(double* x, double* y, double* z, double* bx, double* by, double* bz);
extern "C" void recalc_08_(int* iyr, int* idy, int* ihr, int* imn, int* isc, double* vgsex, double* vgsey, double* vgsez);

// Bmodel:
void bmodel_dipole(double* x_in, double* B_out);
void dipole_geo(int itime_in[2], double x_in[3], double b_out[3]);
void dipole_sm(int itime_in[2], double x_in[3], double b_out[3]);
void bmodel(int itime_in[2], double x_in[3], double tsyg_params[10],
            int use_IGRF, int use_tsyg, int recalc, double b_out[3]);

int trace_fieldline(int itime_in[2], double x_in[3], double x_out[TRACER_MAX][3], 
                    double ds_in, int use_IGRF, int use_tsyg, double tsyg_params[10]);
void init_igrf(int itime_in[2]);
void igrf_geo(double x_in[3], double b_out[3]);
void igrf_mag_cart(int itime_in[2], double x_in[3], double b_out[3], bool recalc);
extern "C" void t04_s_(double* IOPT, double tsyg_params[10], float* PS, 
                      float* X, float* Y, float* Z,
                      float* BX, float* BY, float* BZ);

void load_TS05_params(int itime[2], double TS05_params[10], double VG[3]);

void dump_fieldlines(int itime_in[2], int n_lats, int n_lons, int use_IGRF, int use_tsyg, string filename);


extern "C" {
    extern struct {
         double ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,
                SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33,
                E11,E21,E31,E12,E22,E32,E13,E23,E33;
    } geopack1_;
}

#endif
