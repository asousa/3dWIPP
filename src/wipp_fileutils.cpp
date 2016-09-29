#include <damping.h>

using namespace std;
using namespace Eigen;


map<int, rayF> read_rayfile(string fileName) {
   // Read data from a ray file. Returns a map of structs;
   // 1 key per individual ray in the file.
   //  Each structure contains:

   //  Coordinates:
   //  ray.pos   (n x 3 double)   Position in SM cartesian coordinates
   //  ray.vprel (n x 3 double)   Phase velocity (relative to c)
   //  ray.vgrel (n x 3 double)   Group velocity (relative to c)
   //  ray.n     (n x 3 double)   Index of refraction along each axis
   //  ray.B0    (n x 3 double)   Background magnetic field vector

   //  Constants:
   //  ray.w         double                    Wave angular frequency (radians)
   //  ray.stopcond  double                    Reason for raytracing termination (see docs)
   //  ray.nspec     double                    number of species used in plasmasphere model
    
   //  Plasma parameters:
   //  ray.qs    (nspec x 1 double)            Species charge in coulombs
   //  ray.ms    (nspec x 1 double)            Species mass in kg
   //  ray.Ns    (n x nspec double)            Species number density in m^-3
   //  ray.nus   (n x nspec double)            Species collision frequencies in s^-1

    FILE * filePtr;
    ifstream file;
    string token;
    string line;
    rayF rf;
    int linecounter = 0;
    // int vec_length = 100;

    // Temp doubles for fscanf:
    double ray_num, stopcond, time;
    // double pos_x, pos_y, pos_z;
    // double vprel_x, vprel_y, vprel_z;
    // double vgrel_x, vgrel_y, vgrel_z;
    // double n_x, n_y, n_z;
    // double B0_x, B0_y, B0_z;
    double w, nspec;
    
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

            // // Build an istream that holds the input string
            iss.str(line);

            // // Iterate over the istream, using >> to grab doubles
            // // and push_back to store them in the vector
            copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(v));

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
                Nsv.push_back(v[20 + 2*nspec + i]);
                nuv.push_back(v[20 + 3*nspec + i]);
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

        }  // Parse loop
        cout << "Total lines: " << linecounter << "\n";

    } else {   // Couldn't open the file
        cout << "Couldn't open file\n";
    }

    return raylist;
} // Closing Main.


