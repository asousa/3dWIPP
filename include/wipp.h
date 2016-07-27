// wipp.h
#ifndef wipp_H
#define wipp_H

#include <vector>
#include <map>
using namespace std;

// Structure for holding an entire rayfile (yuge)
typedef struct rayfile_struct {
    int size;
    float w;                        // frequency (angular)
    float stopcond;                   // Stop condition
    float nspec;                      // Number of species in plasmasphere model
    // float ray_num;         // Ray index number
    vector <float> time;            // Group time
    vector <float> pos_x;
    vector <float> pos_y;
    vector <float> pos_z;
    vector <float> vprel_x;
    vector <float> vprel_y;
    vector <float> vprel_z;
    vector <float> vgrel_x;
    vector <float> vgrel_y;
    vector <float> vgrel_z;
    vector <float> n_x;
    vector <float> n_y;
    vector <float> n_z;
    vector <float> B0_x;
    vector <float> B0_y;
    vector <float> B0_z;

    // Variable-length stuff (depending on number of constituents in model)
    vector <vector <float> > qs;    // species charge
    vector <vector <float> > ms;    // species mass
    vector <vector <float> > Ns;    // Uh...
    vector <vector <float> > nus;   // yeah

} rayF;

// rayfile loader
map<int, rayF> read_rayfile(string fileName);

#endif