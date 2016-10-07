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
    ray.n     (n x 3 double)   Index of refraction along each axis
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
    cout << "reading " << fileName << "...\n";

    file.open(fileName.c_str());

    if (file.is_open()) {
        cout << "Successfully opened " << fileName << "\n";
                
        while (getline(file, line)) {
            vector <double> v;
            istringstream iss;

            vector <double> qsv;
            vector <double> msv;
            vector <double> Nsv;
            vector <double> nuv;
            vector <double> B0, n, pos, vprel, vgrel;
            // int a, b, c, d, e, f, g, h;
            // char hh[100];
            // if (line[0] == 'H') {

            //     sscanf(line.c_str(), "%*s %u %u %u %u %u %u %u %u", &ray_num,&iyr,&idoy,&isec,&in_radius,&in_lat,&in_lon, &in_w);
            //     // ray_num = double(a);

            //     printf("Found header for ray %u\n",ray_num);
            //     raylist.insert(make_pair(ray_num, rf));

            //     raylist[ray_num].iyr  = iyr;
            //     raylist[ray_num].idoy = idoy;
            //     raylist[ray_num].isec = isec;
            //     raylist[ray_num].in_radius = in_radius;
            //     raylist[ray_num].in_lat = in_lat;
            //     raylist[ray_num].in_lon = in_lon;
            //     raylist[ray_num].in_w = in_w;

            //     printf("in coords: %u, %u, %u\n",in_radius, in_lat, in_lon);


            //     // printf("read vals: %i %i %i\n",ray_num, iyr, idoy);
            //     // print_vector(v);
            //     // ray_num = v[0];
            //     // iyr  = v[1];
            //     // idoy = v[2];
            //     // isec = v[3];
            //     // in_radius= v[4];
            //     // in_lat = v[5];
            //     // in_lon = v[6];
            //     // in_w   = v[7];

            //     // cout << "ray_num: " << ray_num << "\n";

            // } else {
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
                    printf("Adding new entry for ray number %u\n",ray_num);
                    raylist.insert(make_pair(ray_num, rf));
                    // Set single-value elements
                    raylist[ray_num].w = w;
                    raylist[ray_num].nspec = nspec;
                    raylist[ray_num].stopcond = stopcond;
                }

                // if (raylist.count(ray_num) == 1) {
                //     raylist[ray_num].w = w;
                //     raylist[ray_num].nspec = nspec;
                //     raylist[ray_num].stopcond = stopcond;
                // }

                raylist[ray_num].time.push_back(v[2]);

                pos.push_back(v[3]);
                pos.push_back(v[4]);
                pos.push_back(v[5]);
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
                
                linecounter++;

            // }   // Header check

        }  // Parse loop
    
        cout << "Total lines: " << linecounter << "\n";

    } else {   // Couldn't open the file

        cout << "Something's Fucky\n";
    }

    // return 0; // Return statement.
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