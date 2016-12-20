// #include <consts.h>
#include <wipp.h>

using namespace std;

map <int, rayF> read_rayfile(string fileName)
/* Read data from a ray file. Returns a map of structs;
   1 key per individual ray in the file.
    Each structure contains:

    Coordinates:
    ray.pos   (n x 3 double)   Position in SM cartesian coordinates
    ray.vprel (n x 3 double)   Phase velocity (relative to c)
    ray.vgrel (n x 3 double)   Group velocity (relative to c)
    ray.n     (n x 3 double)   Refractive index vector
    ray.B0    (n x 3 double)   Background magnetic field vector

    Constants:
    ray.w         double                    Wave angular frequency (radians)
    ray.stopcond  double                    Reason for raytracing termination (see docs)
    ray.nspec     double                    number of species used in plasmasphere model
    
    Plasma parameters:
    ray.qs    (nspec x 1 double)            Species charge in coulombs
    ray.ms    (nspec x 1 double)            Species mass in kg
    ray.Ns    (n x nspec double)            Species number density in m^-3
    ray.nus   (n x nspec double)            Species collision frequencies in s^-1

*/

{ 
    FILE * filePtr;
    ifstream file;
    string token;
    string line;
    rayF rf;
    int linecounter = 0;
    // int vec_length = 100;

    // Temp doubles for fscanf:
    double stopcond, time;
    int ray_num;
    // double pos_x, pos_y, pos_z;
    // double vprel_x, vprel_y, vprel_z;
    // double vgrel_x, vgrel_y, vgrel_z;
    // double n_x, n_y, n_z;
    // double B0_x, B0_y, B0_z;
    double w, nspec;
    int iyr, idoy, isec;
    int in_radius, in_lat, in_lon, in_w;
    
    map <int, rayF > raylist;
    // vector <double> v;
    // istringstream iss;

    // vector <double> qsv;
    // vector <double> msv;
    // vector <double> Nsv;
    // vector <double> nuv;


    // Say hi
    if (DEBUG) {cout << "reading " << fileName << "...\n";}

    file.open(fileName.c_str());

    if (file.is_open()) {
        if (DEBUG) {cout << "Successfully opened " << fileName << "\n";}
                
        while (getline(file, line)) {
            vector <double> v;
            istringstream iss;

            vector <double> qsv;
            vector <double> msv;
            vector <double> Nsv;
            vector <double> nuv;
            vector <double> B0, n, pos, vprel, vgrel;

            // // Build an istream that holds the input string
            iss.str(line);

            // // Iterate over the istream, using >> to grab doubles
            // // and push_back to store them in the vector
            copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(v));

            // for (vector<double>::iterator it = v.begin(); it != v.end(); ++it)
            //     printf("%g ",*it);
            // cout << "\n";

            // single-valued parameters
            ray_num = int(v[0]);
            stopcond = v[1];
            w = v[18];
            nspec = v[19];


            // Start a new entry if not in the dictionary already:
            if (raylist.count(ray_num) == 0) {
                // printf("Adding new entry for ray number %u\n",ray_num);
                raylist.insert(make_pair(ray_num, rf));
                // Set single-value elements
                raylist[ray_num].w = w;
                raylist[ray_num].nspec = nspec;
                raylist[ray_num].stopcond = stopcond;
            }

            raylist[ray_num].time.push_back(v[2]);
            // Position (convert to earth radii)
            pos.push_back(v[3]/R_E);
            pos.push_back(v[4]/R_E);
            pos.push_back(v[5]/R_E);
            raylist[ray_num].pos.push_back(pos);

            vprel.push_back(v[6]);
            vprel.push_back(v[7]);
            vprel.push_back(v[8]);
            raylist[ray_num].vprel.push_back(vprel);
            
            vgrel.push_back(v[9]);
            vgrel.push_back(v[10]);
            vgrel.push_back(v[11]);
            raylist[ray_num].vgrel.push_back(vgrel);

            n.push_back(v[12]);
            n.push_back(v[13]);
            n.push_back(v[14]);
            raylist[ray_num].n.push_back(n);

            B0.push_back(v[15]);
            B0.push_back(v[16]);
            B0.push_back(v[17]);
            raylist[ray_num].B0.push_back(B0);

            for (int i = 0; i < nspec; i++) {
                // cout << i << ' ';
                // qsv.push_back(v[20 + 0*nspec + i]);
                // msv.push_back(v[20 + 1*nspec + i]);
                Nsv.push_back(v[20 + 2*nspec + i]);
                nuv.push_back(v[20 + 3*nspec + i]);

                // cout << Nsv[i] << ' ';
            }
            raylist[ray_num].Ns.push_back(Nsv);
            raylist[ray_num].nus.push_back(nuv);

            // Qs, Ms are constants for the whole ray. 
            // I mean why wouldn't they be, right, right
            if (raylist[ray_num].qs.size() == 0) {
                for (int i = 0; i < nspec; i++) {
                    raylist[ray_num].qs.push_back(v[20 + 0*nspec + i]);
                    raylist[ray_num].ms.push_back(v[20 + 1*nspec + i]);
                }
            }

            // Do we have damping data? if so load it too.
            if (v.size() > (20 + 4*nspec) ) {
                // cout <<"damping " << v[20+4*nspec] <<"\n";
                raylist[ray_num].damping.push_back(v[20+4*nspec]);
            }
            
            linecounter++;

        }  // Parse loop
    
        if (DEBUG) {cout << "Total lines: " << linecounter << "\n";}

    } else {   // Couldn't open the file

        cout << "Failed to open file at " << fileName << "\n";
    }


    // // Go thru raylist and set start value, length, Stix parameters, etc
    // for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter){
    
    //     ray = &(iter->second);

    //     start_pos = &(ray->pos[0].data()[0]);

    //     // Get magnetic lat:
    //     sm_to_mag_d_(itime_in, start_pos, tmp_coords2);
    //     cart_to_pol_d_(tmp_coords2, &maglat0, &maglon0, &magrad0);
    //     maglat0 = R2D*maglat0; maglon0 = R2D*maglon0;
    //     printf("MAG lat: %g lon: %g alt: %g\n",maglat0,maglon0,magrad0);
        
    //     // Get power scaling:
    //     // (still need to multiply by space + freq bin sizes)
    //     ray->inp_pwr = input_power_scaling(flash_pos_sm, start_pos, maglat0, iter->second.w, flash_I0);
    //     ray->in_radius = magrad0;
    //     ray->in_lat = maglat0;
    //     ray->in_lon = maglon0;

    //     // Calculate Stix parameters:
    //     calc_stix_parameters(ray);
    // }





    return raylist;
} // Closing Main.

map <int, rayF> read_dampfile(string fileName)
/* Read data from a ray file. Returns a map of structs;
   1 key per individual ray in the file.
    Each structure contains:

    Coordinates:
    ray.pos   (n x 3 double)   Position in SM cartesian coordinates
    ray.vprel (n x 3 double)   Phase velocity (relative to c)
    ray.vgrel (n x 3 double)   Group velocity (relative to c)
    ray.n     (n x 3 double)   Refractive index vector
    ray.B0    (n x 3 double)   Background magnetic field vector

    Constants:
    ray.w         double                    Wave angular frequency (radians)
    ray.stopcond  double                    Reason for raytracing termination (see docs)
    ray.nspec     double                    number of species used in plasmasphere model
    
    Plasma parameters:
    ray.qs    (nspec x 1 double)            Species charge in coulombs
    ray.ms    (nspec x 1 double)            Species mass in kg
    ray.Ns    (n x nspec double)            Species number density in m^-3
    ray.nus   (n x nspec double)            Species collision frequencies in s^-1

*/

{ 
    FILE * filePtr;
    ifstream file;
    string token;
    string line;
    rayF rf;
    int linecounter = 0;
    // int vec_length = 100;

    // Temp doubles for fscanf:
    double stopcond, time;
    int ray_num;
    double w, nspec;
    int iyr, idoy, isec;
    int in_radius, in_lat, in_lon, in_w;
    
    map <int, rayF > raylist;

    // Say hi
    if (DEBUG) {cout << "reading " << fileName << "...\n";}

    file.open(fileName.c_str());

    if (file.is_open()) {
        if (DEBUG) {cout << "Successfully opened " << fileName << "\n";}
                
        while (getline(file, line)) {
            vector <double> v;
            istringstream iss;

            // // Build an istream that holds the input string
            iss.str(line);

            // // Iterate over the istream, using >> to grab doubles
            // // and push_back to store them in the vector
            copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(v));

            ray_num = int(v[0]);

            // Start a new entry if not in the dictionary already:
            if (raylist.count(ray_num) == 0) {
                raylist.insert(make_pair(ray_num, rf));
            }

            raylist[ray_num].time.push_back(v[1]);        
            raylist[ray_num].damping.push_back(v[2]);
                
            linecounter++;

        }  // Parse loop
    
        if (DEBUG)  {cout << "Total lines: " << linecounter << "\n";}

    } else {   // Couldn't open the file

        cout << "Failed to open file at " << fileName << "\n";
    }

    return raylist;
} // Closing Main.





void write_rayfile(string fileName, map <int, rayF> raylist) {
    rayF ray;
    FILE * file;

    cout << "writing " << fileName << "...\n";
    
    file = fopen(fileName.c_str(), "w");    
    if (file != NULL) {
        cout << "Successfully opened " << fileName << "\n";

        // Loop over entries in raylist:
        for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter){
            ray = iter->second;
            cout << ray.time.size() << "\n";

            // Loop over entries in each ray:
            for (int x=0; x < ray.time.size(); x++) {
                // Re-print everything in the same format as read in:
                // (Originally specified in the fortran raytracer)
                fprintf(file, "%d  %g  %.15e  ", iter->first, ray.stopcond, ray.time[x]);
                fprintf(file, "%.15e  %.15e  %.15e  ", ray.pos[x][0], ray.pos[x][1], ray.pos[x][2]);
                fprintf(file, "%.15e  %.15e  %.15e  ",ray.vprel[x][0], ray.vprel[x][1], ray.vprel[x][2]);
                fprintf(file, "%.15e  %.15e  %.15e  ",ray.vgrel[x][0], ray.vgrel[x][1], ray.vgrel[x][2]);
                fprintf(file, "%.15e  %.15e  %.15e  ",ray.n[x][0],     ray.n[x][1],     ray.n[x][2]);
                fprintf(file, "%.15e  %.15e  %.15e  ",ray.B0[x][0],    ray.B0[x][1],    ray.B0[x][2]);
                fprintf(file, "%.15e  %g  ",ray.w, ray.nspec);
                for (int i = 0; i < ray.nspec; i++) {
                    fprintf(file, "%.15e  ",ray.qs[i]);
                }
                for (int i = 0; i < ray.nspec; i++) {
                    fprintf(file, "%.15e  ",ray.ms[i]);
                }
                for (int i = 0; i < ray.nspec; i++) {
                    fprintf(file, "%.15e  ",ray.Ns[x][i]);
                }
                for (int i = 0; i < ray.nspec; i++) {
                    fprintf(file, "%.15e  ",ray.nus[x][i]);
                }

                // If we have damping calculations, let's append them too:
                if (ray.damping.size() > 0) {

                    fprintf(file, "%.15e  ",ray.damping[x]);
                }

                fprintf(file, "\n");

            } // timestep
        } // raylist

    } else {   // Couldn't open the file
        cout << "Something's Fucky\n";
    }   

}



void load_TS05_params(int itime_in[2], double TS_params[10], double VG[3]) {

    ostringstream inp_filename;
    ifstream file;
    vector <double> v;
    vector <double> v_prev;
    istringstream iss;
    string line;

    double total_mins = 0;
    double dtL, dtR, dt;

    inp_filename << "data/TS04/";

    int yr = (int)(itime_in[0]/1000);
    int doy  = itime_in[0] - 1000*yr;
    int isec = itime_in[1]*1e-3; 
    int hr = (int)(isec/3600);
    int mn = (int)((isec - hr*3600)/60);

    int target_total_mins = mn + hr*60 + doy*24*60;

    // cout << "target total: " << target_total_mins << "\n";
    // cout << "yr: " << yr << " doy: " << doy << " hr: " << hr << " min: " << mn << "\n";
    if (1995 <= yr <=2015) {
        inp_filename << yr << "_OMNI_5m_with_TS05_variables.dat";
        // cout << "inp filename: " << inp_filename.str() << "\n";

        file.open(inp_filename.str().c_str());
        if (file.is_open()) {
            while (getline(file, line)) {

                v_prev = v;
                // v = {0};
                v.clear();
                iss.clear();
                iss.str(line);

                 // // Iterate over the istream, using >> to grab doubles
                // // and push_back to store them in the vector
                copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(v));
                // cout << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << " " << v[4] << "\n";
                total_mins = v[3] + v[2]*60 + v[1]*24*60;
                // cout << "targ mins: " << target_total_mins <<  " total mins: " << total_mins << "\n";
               if (total_mins >= target_total_mins)  {
                    // Interpolate between this row and previous row
                    if (v_prev.size() > 0) {

                        // print_vector(v);
                        // seconds from each segment
                        dtL = target_total_mins - v_prev[3] - v_prev[2]*60 - v_prev[1]*24*60;
                        dtR = total_mins - target_total_mins;
                        dt  = dtR + dtL;
                        // cout << "L: " << dtL << " R: " << dtR << " Tot: " << dt << "\n";
                        // cout << "Tilt: " << v[15] << " " << v_prev[15] << "\n";

                        TS_params[0] = (dtL/dt)*v[16]  +   (dtR/dt)*v_prev[16];   // Pdyn
                        TS_params[1] = 0               +   0;                     // Dst
                        TS_params[2] = (dtL/dt)*v[5]   +   (dtR/dt)*v_prev[5];    // ByIMF
                        TS_params[3] = (dtL/dt)*v[6]   +   (dtR/dt)*v_prev[6];    // BzIMF
                        TS_params[4] = (dtL/dt)*v[17]  +   (dtR/dt)*v_prev[17];   // W1
                        TS_params[5] = (dtL/dt)*v[18]  +   (dtR/dt)*v_prev[18];   // W2
                        TS_params[6] = (dtL/dt)*v[19]  +   (dtR/dt)*v_prev[19];   // W3
                        TS_params[7] = (dtL/dt)*v[20]  +   (dtR/dt)*v_prev[20];   // W4
                        TS_params[8] = (dtL/dt)*v[21]  +   (dtR/dt)*v_prev[21];   // W5
                        TS_params[9] = (dtL/dt)*v[22]  +   (dtR/dt)*v_prev[22];   // W6

                        VG[0]        = (dtL/dt)*v[7]   +   (dtR/dt)*v_prev[7];
                        VG[1]        = (dtL/dt)*v[8]   +   (dtR/dt)*v_prev[8];
                        VG[2]        = (dtL/dt)*v[9]   +   (dtR/dt)*v_prev[9];
                        // TS_params order:
                        // [Pdyn Dst ByIMF BzIMF W1 W2 W3 W4 W5 W6]

                    } else {
                        // first row.
                        TS_params[0] = v[16];   // Pdyn
                        TS_params[1] = 0;       // Dst
                        TS_params[2] = v[5];    // ByIMF
                        TS_params[3] = v[6];    // BzIMF
                        TS_params[4] = v[17];   // W1
                        TS_params[5] = v[18];   // W2
                        TS_params[6] = v[19];   // W3
                        TS_params[7] = v[20];   // W4
                        TS_params[8] = v[21];   // W5
                        TS_params[9] = v[22];   // W6

                    }
                    break;

                } 
            }
            file.close();
        }


    } else {
        cout << "Out of bounds of TS files!\n";
    }
}

// This version for ascii files:
void write_p_array(double arr[NUM_E][NUM_TIMES], string filename) {
    FILE* file;
    
    int writecount;
    
    cout << "file shape: " << NUM_E << ", " << NUM_TIMES << "\n";
    file = fopen(filename.c_str(),"wb");

    cout << "sizeof (arr): " << sizeof(arr) << "\n";
    if (file == NULL) {
        cout << "Failed to open file " << filename << "\n";
    } else {

        // // ASCII-formatted:
        // for (int x=0; x < NUM_E; x++) {
        //     for (int y=0; y < NUM_TIMES; y++) {
        //        fprintf(file, "%g ",arr[x][y]);
        //     }    
        // }

        // // Binary-formatted:
        // // setbuf(file, NULL);
        // for(int ei=0; ei<NUM_E; ei++) {
        //     for(int ti=0; ti<NUM_TIMES; ti++) {
        //         fwrite(&arr[ei][ti], sizeof(float), 1, file);
        //     }
        // }
        fwrite(arr, NUM_E*NUM_TIMES*sizeof(double), 1, file);

    }
    fclose(file);
}



void read_p_array(double arr[NUM_E][NUM_TIMES], string filename) {
    ifstream file;
    vector <double> v;
    istringstream iss;
    string line;

    file.open(filename.c_str());
    if (file.is_open()) {
        while (getline(file, line)) {

            v.clear();
            iss.clear();
            iss.str(line);

            copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(v));

            if (v.size() != NUM_E*NUM_TIMES) {
                cout << "size mismatch!\n";
            } else {
                for (int row=0; row < NUM_E; ++row) {
                    for (int col=0; col < NUM_TIMES; ++col){
                        arr[row][col] = v[row*NUM_TIMES + col];
                    }
                }
            }
        }
    }
}




vector<cellT> load_crossings(int itime_in[2], string filename) {
    ifstream file;
    istringstream iss;
    string line;
    double x_in[3];
    double x_out[3], x_sm[3];

    vector<cellT> outs;

    double Lsh, lat, t, f, pwr, psi, mu, stixP, stixR, stixL;
    cout << "Fuckin jesus \n";
    file.open(filename.c_str());
    if (file.is_open()) {
        while (getline(file, line)) {
            vector <double> v;
            istringstream iss;

            // cout << line;
            cellT ray = {};


            iss.str(line);
            
            //  // // Iterate over the istream, using >> to grab doubles
            // // // and push_back to store them in the vector
            copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(v));
            
            // print_vector(v);

            Lsh = v[0]; lat = v[1]; t = v[2]; f = v[3]; pwr = v[4]; psi = v[5]; mu = v[6];
            stixP = v[7]; stixR = v[8]; stixL = v[9];
            // Get coordinates in SM:
            double radius = Lsh*pow(cos(D2R*lat),2);  
            x_in = {radius, D2R*lat, 0};
            sphcar(x_in, x_out);
            mag_to_sm_d_(itime_in, x_out, x_sm);


            ray.pos = Map<Vector3d>(x_sm, 3,1);
            ray.pwr = pwr;
            ray.mu  = mu;
            ray.t   = t;
            ray.f   = f;
            ray.psi = psi;
            ray.stixP = stixP;
            ray.stixR = stixR;
            ray.stixL = stixL;
            ray.num_rays = 1;
            ray.Lsh = Lsh;
            ray.lat = lat;

            outs.push_back(ray);
        }
    }
    return outs;
}


void get_available_rays(string raypath, vector <vector<double> > *data) {
// Recursively search through a file tree, and pull out all frequencies, latitudes, and longitudes.
    ostringstream curpath;
    float freq, lat, lon;
    DIR *dp;
    struct dirent *ep;    
    cout << raypath << "\n";
    dp = opendir (raypath.c_str());
    if (dp != NULL) {
        while (ep = readdir(dp)) {
            if (strcmp(ep->d_name,".") && strcmp(ep->d_name, "..")) {   /* ignore up directories */
                curpath.str(""); curpath.clear();
                curpath << raypath << "/" <<ep->d_name;
                // Try to recurse
                get_available_rays(curpath.str(), data);
            }
            if (sscanf(ep->d_name, "ray_%g_%g_%g.ray",&freq, &lat, &lon)) {
                vector<double> row;
                row.push_back(freq); row.push_back(lat); row.push_back(lon);
                data->push_back(row);
            }
        }
    }
}




int check_memory_usage() {
    int tSize = 0, resident = 0, share = 0;
    ifstream buffer("/proc/self/statm");
    buffer >> tSize >> resident >> share;
    buffer.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    double rss = resident * page_size_kb;
    cout << "RSS - " << rss << " kB\n";

    double shared_mem = share * page_size_kb;
    cout << "Shared Memory - " << shared_mem << " kB\n";

    cout << "Private Memory - " << rss - shared_mem << "kB\n";
    return 0;
}


// static int ftw_callback(const char *fpath, const struct stat *sb, int typeflag) {
//     float freq, lat, lon;
//     const char* ms;
//     if (typeflag == FTW_F) {

//         ms = strstr(fpath, "ray_");
//         if (ms != NULL) {
//             // cout << fpath <<"\n";
//             // cout << ms ;
//             sscanf(ms, "ray_%g_%g_%g.ray", &freq, &lat, &lon);

        
//         // // sscanf(fpath, "%sray_%g_%g_%g.ray",&str, &freq, &lat, &lon);
//             cout << " " << freq << " " << lat << " " << lon << "\n";
// }        
//         // cout << fpath << "\n";
//     }
//     return 0;
// }