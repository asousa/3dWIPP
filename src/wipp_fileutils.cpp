#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
// #include <string>
#include <vector>
#include <consts.h>
#include <wipp.h>


// #include <boost/algorithm/string.hpp>

using namespace std;


// Structure for holding an entire rayfile (yuge)
typedef struct {
    int size;
    vector <float> ray_num;         // Ray index number
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
} rayF;

// Prototypes:
// void resize_rayfile(rayF &ray, int newsize);


int read_rayfile(string fileName) 
{ 
    FILE * filePtr;
    ifstream file;
    string token;
    string line;

    rayF rf;
    int linecounter = 0;
    int vec_length = 100;

    // Temp floats for fscanf:
    float ray_num;
    float time;
    float pos_x, pos_y, pos_z;
    float vprel_x, vprel_y, vprel_z;
    float vgrel_x, vgrel_y, vgrel_z;
    float n_x, n_y, n_z;
    float B0_x, B0_y, B0_z;

    float tmpfloat;
    
    vector <rayF> raylist;

    // Say hi
    cout << "reading " << fileName << "...\n";

    file.open(fileName.c_str());

    if (file.is_open()) {
        cout << "Successfully opened " << fileName << "\n";
                
        while (getline(file, line)) {

            sscanf(line.c_str(), "%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e",
                &ray_num, &time, 
                &pos_x,   &pos_y,   &pos_z,
                &vprel_x, &vprel_y, &vprel_z,
                &vgrel_x, &vgrel_y, &vgrel_z,
                &n_x,     &n_y,     &n_z,
                &B0_x,    &B0_y,    &B0_z);


            if (raylist.size() <= ray_num) {
                cout << "Adding to raylist\n";
                cout << ray_num;
                raylist.resize(ray_num + 1);
            }

            // Push temps to vectors
            raylist[ray_num].ray_num.push_back(ray_num);
            raylist[ray_num].time.push_back(time);

            raylist[ray_num].pos_x.push_back(pos_x);
            raylist[ray_num].pos_y.push_back(pos_y);
            raylist[ray_num].pos_z.push_back(pos_z);

            raylist[ray_num].vprel_x.push_back(vprel_x);
            raylist[ray_num].vprel_y.push_back(vprel_y);
            raylist[ray_num].vprel_z.push_back(vprel_z);
            
            raylist[ray_num].vgrel_x.push_back(vgrel_x);
            raylist[ray_num].vgrel_y.push_back(vgrel_y);
            raylist[ray_num].vgrel_z.push_back(vgrel_z);

            raylist[ray_num].n_x.push_back(n_x);
            raylist[ray_num].n_y.push_back(n_y);
            raylist[ray_num].n_z.push_back(n_z);
            
            raylist[ray_num].B0_x.push_back(B0_x);
            raylist[ray_num].B0_y.push_back(B0_y);
            raylist[ray_num].B0_z.push_back(B0_z);
            
            // cout << rf.ray_num[linecounter] << ' ';
            
            linecounter++;

        }  // Parse loop
    
        cout << "Total lines: " << linecounter << "\n";

        // for (int i=0; i < raylist.size(); i++) {
        //     for (int j=0; j < raylist[i].ray_num.size(); j++) {

        //         cout << raylist[i].ray_num[j] << " ";
        //     }
        //     cout << "\n";
        // }

    } else {   // Couldn't open the file

        cout << "Something's Fucky\n";
    }

    return 0; // Return statement.
} // Closing Main.







