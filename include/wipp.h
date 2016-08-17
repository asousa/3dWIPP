// wipp.h
#ifndef wipp_H
#define wipp_H

#include <vector>
#include <map>

using namespace std;

// Structure for holding an entire rayfile (yuge)
typedef struct rayfile_struct {
    int size;
    double w;                        // frequency (angular)
    double stopcond;                   // Stop condition
    double nspec;                      // Number of species in plasmasphere model
    // double ray_num;         // Ray index number
    vector <double> time;            // Group time
    // vector <double> pos_x;
    // vector <double> pos_y;
    // vector <double> pos_z;
    // vector <double> vprel_x;
    // vector <double> vprel_y;
    // vector <double> vprel_z;
    // vector <double> vgrel_x;
    // vector <double> vgrel_y;
    // vector <double> vgrel_z;
    // vector <double> n_x;
    // vector <double> n_y;
    // vector <double> n_z;
    vector <vector <double> > pos;
    vector <vector <double> > vprel;
    vector <vector <double> > vgrel;
    vector <vector <double> > n;
    vector <vector <double> > B0;


    // vector <double> B0_x;
    // vector <double> B0_y;
    // vector <double> B0_z;

    // Variable-length stuff (depending on number of constituents in model)
    vector <double> qs;    // species charge
    vector <double> ms;    // species mass
    vector <vector <double> > Ns;    // number density of species (m^-3)
    vector <vector <double> > nus;   // collision frequencies

} rayF;

// rayfile loader
map<int, rayF> read_rayfile(string fileName);

// Landau damping
void damping_ngo(rayF rayfile);

// Math functions
double l2_norm(vector<double> u);
vector<double> scalar_multiply(vector<double> u, double v);
double dot_product(vector<double>u, vector<double>v);
vector<double> add(vector<double>u, vector<double> v);

// Helpers
void print_vector(vector<double> u);

#endif