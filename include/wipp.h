// wipp.h
#ifndef wipp_H
#define wipp_H

#include <vector>
#include <map>

#include <mat.h>

using namespace std;

// Structure for holding an entire rayfile (yuge)
typedef struct rayfile_struct {
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

    // Variable-length stuff (depending on number of constituents in model)
    vector <double> qs;    // species charge
    vector <double> ms;    // species mass
    vector <vector <double> > Ns;    // number density of species (m^-3)
    vector <vector <double> > nus;   // collision frequencies
    vector <double> damping;
} rayF;

// rayfile loader
map<int, rayF> read_rayfile(string fileName);

// Landau damping
void damping_ngo(rayF &rayfile);

// Math functions
double l2_norm(vector<double> u);
vector<double> scalar_multiply(vector<double> u, double v);
double dot_product(vector<double>u, vector<double>v);
vector<double> add(vector<double>u, vector<double> v);

// Helpers
void print_vector(vector<double> u);


// Porting the Damping Code:
void damping_foust(rayF &rayfile);

double kp_to_pp(double kp);
void polyfit(const vector<double> &xv, const vector<double> &yv, vector<double> &coeff, int order);

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
