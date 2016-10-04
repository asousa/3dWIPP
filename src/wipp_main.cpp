#include <wipp.h>

using namespace std;
using namespace Eigen;

// External functions we'll use (libxformd for coordinate transforms)
// extern "C" void sm_to_geo_d_(int* itime, double* x_in, double* x_out);
// extern "C" void sm_to_mag_d_(int* itime, double* x_in, double* x_out);
// extern "C" void geo_to_sm_d_(int* itime, double* x_in, double* x_out);
// extern "C" void cart_to_pol_d_(double* x_in, double* lat, double* lon, double* radius);
// extern "C" void pol_to_cart_d_(double* lat, double* lon, double* radius, double* x_out);

// lib_onera_desp (the Cospar IRBEM library, which is the internal support for SpacePy):
extern "C" void sm2geo1_(int* iyr,int* idoy, double* secs, double* xSM, double* xGEO);
extern "C" void geo2sm1_(int* iyr,int* idoy, double* secs, double* xGEO, double* xSM);

// cartesian - spherical (trig terms in degrees!)
extern "C" void car_sph_(double* xCAR, double* r, double* lat, double* loni);
extern "C" void sph_car_(double* r, double* lat, double* loni, double* xCAR);
 
int main(int argc, char *argv[]) 
{
    map <int, rayF> raylist;
    map <int, VectorXd> damplist;


    double x_in[2];
    double x_out[2];

    string inpFileName;
    string outFileName;
    
    int itime_in[2];
    
    string itime_str;

    double tmp_coords[2];
    double tmp_in[2];
    double lat0, lon0, rad0;

    FILE * outputFile;

    // Default parameters:
    inpFileName = "input.ray";
    outFileName = "output.ray";

    int iyr = 2010;     // Year
    int idoy= 155;      // Day of year
    double isec= 1.0*60*60*12; // Seconds into day (UTC)

     // Parse input arguments:
    int opt = 0;
    while ((opt = getopt(argc, argv, "i:o:t:u:v:a:b:c:")) != -1) {
        switch(opt) {
            case 'i':
            // input filename:
            // cout << "inp filename!\n";
                inpFileName = (string) optarg;
                break;
            case 'o':
                outFileName = (string) optarg;
                break;
            case 't':
                iyr = atoi(optarg);
                break;
            case 'u':
                idoy= atoi(optarg);
                break;
            case 'v':
                isec= strtod(optarg, NULL);
                break;
            case 'a':
                x_in[0] = strtod(optarg, NULL);
                break;
            case 'b':
                x_in[1] = strtod(optarg, NULL);
                break;
            case 'c':
                x_in[2] = strtod(optarg, NULL);
                break;
            case '?':
                 printf("\nUnknown option: %s\n",opt);
            break;
        }
    }

    cout << "---- Input Parameters ----\n";
    cout << "input file: " << inpFileName << "\n";
    cout << "output file: " << outFileName << "\n";
    cout << "year: " << iyr << "\n";
    cout << "day: " << idoy << "\n";
    cout << "sec of day: " << isec << "\n";
    cout << "---- 3D WIPP ----\n";

    // tmp_in[0] =-3653337.8925629;
    // tmp_in[1] =-2740585.99220064;
    // tmp_in[2] =4452100.02831447;
    
          
    cout << "x_in: " << x_in[0] << ", " << x_in[1] << ", " << x_in[2] << "\n";

    // sm_to_geo_d_(itime_in, x_in, tmp_coords);
    // cart_to_pol_d_(tmp_coords, &lat0, &lon0, &rad0);
    // printf("ray origin (xformd): %g, %g, %g\n",lat0*R2D, lon0*R2D, rad0/1000.);

    // int iyr = int(floor(itime_in[0]/1000.0));
    // int idoy= int(fmod(itime_in[0], float(iyr)));
    // double secs = itime_in[1]*1e-3;

    cout << "iyr: " << iyr << ", idoy: " << idoy << "\n";
    sm2geo1_(&iyr, &idoy, &isec, x_in, x_out);
    car_sph_(x_out, &rad0, &lat0, &lon0);
    // cart_to_pol_d_(x_out, &lat0, &lon0, &rad0);
    printf("ray origin (other thing): %g, %g, %g\n",lat0, lon0, rad0/1000.);




    // // Load the rayfile:
    // raylist = read_rayfile(inpFileName);

    // vector <vector <double> > start_pos;

    // // Parse rays, make a nice list of the start coords and 
    // for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter){
    //     // printf("damping ray # %d\n",iter->first);
    //     start_pos.push_back(iter->second.pos[0]);
    //     // printf("ray origin: %g, %g, %g\n",iter->second.pos[0][0], iter->second.pos[0][1], iter->second.pos[0][2]);
    //     sm_to_geo_d_(itime_in, &(iter->second.pos[0][0]), tmp_coords);
    //     cart_to_pol_d_(tmp_coords, &lat0, &lon0, &rad0);
    //     printf("lat: %g lon: %g alt: %g\n",lat0*R2D,lon0*R2D,rad0/1000.);
    //     // cout << iter->second.pos[0][0] << ", " << iter->second.pos[0][1] << iter->second.pos[0][2] << "\n";
    // }



    return 0; // Return statement.
} // Closing Main.