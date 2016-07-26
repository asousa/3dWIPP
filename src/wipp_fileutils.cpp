#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
// #include <string>
#include <vector>
#include <map>
#include <consts.h>
#include <wipp.h>


// #include <boost/algorithm/string.hpp>

using namespace std;


// Structure for holding an entire rayfile (yuge)
typedef struct {
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
    // int vec_length = 100;

    // Temp floats for fscanf:
    float ray_num, stopcond, time;
    // float pos_x, pos_y, pos_z;
    // float vprel_x, vprel_y, vprel_z;
    // float vgrel_x, vgrel_y, vgrel_z;
    // float n_x, n_y, n_z;
    // float B0_x, B0_y, B0_z;
    float w, nspec;
    
    map <int, rayF > raylist;
    // vector <float> v;
    // istringstream iss;

    // vector <float> qsv;
    // vector <float> msv;
    // vector <float> Nsv;
    // vector <float> nuv;



    // Say hi
    cout << "reading " << fileName << "...\n";

    file.open(fileName.c_str());

    if (file.is_open()) {
        cout << "Successfully opened " << fileName << "\n";
                
        while (getline(file, line)) {
            vector <float> v;
            istringstream iss;

            vector <float> qsv;
            vector <float> msv;
            vector <float> Nsv;
            vector <float> nuv;

            // v.clear();
            // iss.clear();
            // qsv.clear();
            // msv.clear();
            // Nsv.clear();
            // nuv.clear();

            // --------------- SCANF way -----------------
            // sscanf(line.c_str(), "%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e",
            //     &ray_num, &stopcond, &time, 
            //     &pos_x,   &pos_y,    &pos_z,
            //     &vprel_x, &vprel_y,  &vprel_z,
            //     &vgrel_x, &vgrel_y,  &vgrel_z,
            //     &n_x,     &n_y,      &n_z,
            //     &B0_x,    &B0_y,     &B0_z,
            // //     &w,       &nspec);

            // raylist[ray_num].time.push_back(time);

            // raylist[ray_num].pos_x.push_back(pos_x);
            // raylist[ray_num].pos_y.push_back(pos_y);
            // raylist[ray_num].pos_z.push_back(pos_z);

            // raylist[ray_num].vprel_x.push_back(vprel_x);
            // raylist[ray_num].vprel_y.push_back(vprel_y);
            // raylist[ray_num].vprel_z.push_back(vprel_z);
            
            // raylist[ray_num].vgrel_x.push_back(vgrel_x);
            // raylist[ray_num].vgrel_y.push_back(vgrel_y);
            // raylist[ray_num].vgrel_z.push_back(vgrel_z);

            // raylist[ray_num].n_x.push_back(n_x);
            // raylist[ray_num].n_y.push_back(n_y);
            // raylist[ray_num].n_z.push_back(n_z);
            
            // raylist[ray_num].B0_x.push_back(B0_x);
            // raylist[ray_num].B0_y.push_back(B0_y);
            // raylist[ray_num].B0_z.push_back(B0_z);
            // ------------------------------------------


            // // Try it with stringstream too
            // // Build an istream that holds the input string
            iss.str(line);

            // // Iterate over the istream, using >> to grab floats
            // // and push_back to store them in the vector
            copy(istream_iterator<float>(iss), istream_iterator<float>(), back_inserter(v));

            // for (vector<float>::iterator it = v.begin(); it != v.end(); ++it)
            //     printf("%g ",*it);
            // cout << "\n";

            // single-valued parameters
            ray_num = v[0];
            stopcond = v[1];
            w = v[18];
            nspec = v[19];


            // Start a new entry if not in the dictionary already:
            if (raylist.count(ray_num) == 0) {
                printf("Adding new entry for ray number %g\n",ray_num);

                raylist.insert(make_pair(ray_num, rf));
                // Set single-value elements
                raylist[ray_num].w = w;
                raylist[ray_num].nspec = nspec;
                raylist[ray_num].stopcond = stopcond;
            }

            raylist[ray_num].time.push_back(v[2]);

            raylist[ray_num].pos_x.push_back(v[3]);
            raylist[ray_num].pos_y.push_back(v[4]);
            raylist[ray_num].pos_z.push_back(v[5]);

            raylist[ray_num].vprel_x.push_back(v[6]);
            raylist[ray_num].vprel_y.push_back(v[7]);
            raylist[ray_num].vprel_z.push_back(v[8]);
            
            raylist[ray_num].vgrel_x.push_back(v[9]);
            raylist[ray_num].vgrel_y.push_back(v[10]);
            raylist[ray_num].vgrel_z.push_back(v[11]);

            raylist[ray_num].n_x.push_back(v[12]);
            raylist[ray_num].n_y.push_back(v[13]);
            raylist[ray_num].n_z.push_back(v[14]);
            
            raylist[ray_num].B0_x.push_back(v[15]);
            raylist[ray_num].B0_y.push_back(v[16]);
            raylist[ray_num].B0_z.push_back(v[17]);

            // vector <float> qsv;
            // vector <float> msv;
            // vector <float> Nsv;
            // vector <float> nuv;
            for (int i = 0; i < nspec; i++) {
                // cout << i << ' ';
                qsv.push_back(v[20 + 0*nspec + i]);
                msv.push_back(v[20 + 1*nspec + i]);
                Nsv.push_back(v[20 + 2*nspec + i]);
                nuv.push_back(v[20 + 3*nspec + i]);

                // cout << Nsv[i] << ' ';
            }
            raylist[ray_num].qs.push_back(qsv);
            raylist[ray_num].ms.push_back(msv);
            raylist[ray_num].Ns.push_back(Nsv);
            raylist[ray_num].nus.push_back(nuv);
            
            linecounter++;

        }  // Parse loop
    
        cout << "Total lines: " << linecounter << "\n";

        // Confirm we loaded everything
        for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter)
            {
                float key = iter->first;
                rayF val = iter->second;
                printf("Ray number %g: freq: %g nspec: %g\n",key, val.w, val.nspec);

                // Print out all elements in a vector
                vector<float> vec = val.B0_x;
                for (vector<float>::iterator it = vec.begin(); it != vec.end(); ++it)
                    printf("%g ",*it);
                cout << "\n";

            }

    } else {   // Couldn't open the file

        cout << "Something's Fucky\n";
    }

    return 0; // Return statement.
} // Closing Main.







