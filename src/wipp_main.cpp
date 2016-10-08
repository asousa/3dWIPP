#include <wipp.h>

using namespace std;
using namespace Eigen;

// External functions we'll use (libxformd for coordinate transforms)
// extern "C" void sm_to_geo_d_(int* itime, double* x_in, double* x_out);
// extern "C" void sm_to_mag_d_(int* itime, double* x_in, double* x_out);
// extern "C" void geo_to_sm_d_(int* itime, double* x_in, double* x_out);
// extern "C" void cart_to_pol_d_(double* x_in, double* lat, double* lon, double* radius);
// extern "C" void pol_to_cart_d_(double* lat, double* lon, double* radius, double* x_out);

// // lib_onera_desp (the Cospar IRBEM library, which is the internal support for SpacePy):
// extern "C" void sm2geo1_(int* iyr,int* idoy, double* secs, double* xSM, double* xGEO);
// extern "C" void geo2sm1_(int* iyr,int* idoy, double* secs, double* xGEO, double* xSM);

// // cartesian - spherical (trig terms in degrees!)
// extern "C" void car_sph_(double* xCAR, double* r, double* lat, double* loni);
// extern "C" void sph_car_(double* r, double* lat, double* loni, double* xCAR);
 
int main(int argc, char *argv[]) 
{
    map <int, rayF> raylist;
    map <int, VectorXd> damplist;

    double flash_pos[3];  // Flash position, geographic spherical coords (r, lat, lon)
    double flash_pos_sm[3];
    double x_in[3];
    double x_out[3];

    string inpFileName;
    string outFileName;

    int itime_in[2];
    
    string itime_str;

    double tmp_coords[3];
    double tmp_coords2[3];

    double x[3];
    double lat0, lon0, rad0;
    double maglat0, maglon0, magrad0;

    double dd;

    // Inputs to the IRBEM coordinate transform library
    long ntimes, c1, c2;

    FILE * outputFile;

    // Default parameters:
    inpFileName = "input.ray";
    outFileName = "output.ray";

    long iyr;     // Year
    long idoy;      // Day of year
    double isec; // Seconds into day (UTC)

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
                // isec= strtod(optarg, NULL);
                isec= atoi(optarg);

                break;
            case 'a':
                flash_pos[0] = strtod(optarg, NULL);
                break;
            case 'b':
                flash_pos[1] = strtod(optarg, NULL);
                break;
            case 'c':
                flash_pos[2] = strtod(optarg, NULL);
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



    cout << "flash geo: ";
    print_vector(vector<double>(flash_pos, flash_pos + 3));
    double flash_I0 = 100e3;
    

    int yearday = iyr*1000 + idoy;


    itime_in[0] = yearday;
    itime_in[1] = isec*1e3;

    lat0 = D2R*flash_pos[1];
    lon0 = D2R*flash_pos[2];
    rad0 = flash_pos[0];

    // Get flash position in SM coords:
    pol_to_cart_d_(&lat0, &lon0, &rad0, tmp_coords);
    mag_to_sm_d_(itime_in, tmp_coords, flash_pos_sm);

    cout << "SM (libxformd): ";
    print_vector(vector<double>(flash_pos_sm, flash_pos_sm + 3));


    // Load the rayfile:
    raylist = read_rayfile(inpFileName);

    // vector <double> start_pos;
    double* start_pos;

    rayF* ray;
    // double start_pos_geo[3];

    // Array of pointers to current rays of interest
    // (does this copy or just point? I hope it points.)
    rayF* cur_rays[8];

    // Parse rays, make a nice list of the start coords and frequencies
    for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter){
    
        cout << "Ray origin (SM):\n";

        ray = &(iter->second);

        start_pos = &(ray->pos[0].data()[0]);
        print_vector(ray->pos[0]);


        // Get magnetic lat:
        sm_to_mag_d_(itime_in, start_pos, tmp_coords2);
        cart_to_pol_d_(tmp_coords2, &maglat0, &maglon0, &magrad0);
        maglat0 = R2D*maglat0; maglon0 = R2D*maglon0;
        printf("MAG lat: %g lon: %g alt: %g\n",maglat0,maglon0,magrad0);
        
        // Get power scaling:
        // (still need to multiply by space + freq bin sizes)
        ray->inp_pwr = input_power_scaling(flash_pos_sm, start_pos, maglat0, iter->second.w, flash_I0);
        cout << "dd: " << ray->inp_pwr << "\n";
    }



    // cout << raylist.at(1).inp_pwr << "\n";

    cur_rays[0] = &(raylist.at(1));
    cur_rays[1] = &(raylist.at(2));
    cur_rays[2] = &(raylist.at(3));
    cur_rays[3] = &(raylist.at(4));
    cur_rays[4] = &(raylist.at(5));
    cur_rays[5] = &(raylist.at(6));
    cur_rays[6] = &(raylist.at(7));
    cur_rays[7] = &(raylist.at(8));


    cout << "sanity check: " << cur_rays[7]->inp_pwr << "\n";

    double pos_interp[3];
    double weight_ind[3];



    for (double ii=0; ii <= 1; ii+=0.5) {
        for (double jj=0; jj <= 1; jj+= 0.5) {
            for (double kk=0; kk <= 1; kk+= 0.5) {

                pos_interp = {0.,0.,0.};
                
                // weight_ind = {ii, jj, kk};
                int t_ind = 0;

                interp_vector(cur_rays, ii, jj, kk, t_ind, pos_interp);

                // Get magnetic lat:
                sm_to_mag_d_(itime_in, pos_interp, tmp_coords2);
                cart_to_pol_d_(tmp_coords2, &maglat0, &maglon0, &magrad0);
                maglat0 = R2D*maglat0; maglon0 = R2D*maglon0;
                printf("(%g, %g, %g)\tmag lat, lon: (%g, %g)\n",ii, jj, kk, maglat0,maglon0);

            }
        }
    }
    // for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter){

    //     cout << "dd: " << iter->second.inp_pwr << "\n";
    // }
    return 0; // Return statement.
} // Closing Main.