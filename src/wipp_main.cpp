#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>
#include <iostream>

#include <consts.h>
#include <wipp.h>

// #include <lib/libxformd.a>

// #include <boost/program_options.hpp>
// #include <cassert>

using namespace std;

// External functions we'll use (libxformd for coordinate transforms)
extern "C" void sm_to_geo_d_(int* itime, double* x_in, double* x_out);  
 
int main(int argc, char *argv[]) 
{
    map <int, rayF> raylist;

    double x_in[3];
    double x_out[3];
    int itime;

    // char *fileName;
    string fileName;

    if (argc != 2) {
        fileName = "python/four_adjacent.ray";
    } else {
        fileName = argv[1];
    }




    // fileName= "rayfile.ray";
    cout << "---- 3D WIPP ----\n";
    // Load the rayfile:
    raylist = read_rayfile(fileName);

    itime = 0;
    x_in[0] = raylist[1].pos_x[0];
    x_in[1] = raylist[1].pos_x[1];
    x_in[2] = raylist[1].pos_x[2];

    printf("x_in is: %g, %g, %g\n",x_in[0], x_in[1], x_in[2]);
    sm_to_geo_d_(&itime, x_in, x_out);
    printf("x_out is: %g, %g, %g\n",x_out[0], x_out[1], x_out[2]);


    // // Confirm we loaded everything
    // for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter)
    //     {
    //         float key = iter->first;
    //         rayF val = iter->second;
    //         printf("Ray number %g: freq: %g nspec: %g\n",key, val.w, val.nspec);

    //         // Print out all elements in a vector
    //         vector<float> vec = val.time;
    //         for (vector<float>::iterator it = vec.begin(); it != vec.end(); ++it)
    //             printf("%g ",*it);
    //         cout << "\n";

    //     }


    return 0; // Return statement.
} // Closing Main.