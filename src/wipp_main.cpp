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

    string dumpFileName;

    int itime_in[2];
    
    string itime_str;

    double tmp_coords[3];
    double tmp_coords2[3];

    double x[3];
    double lat0, lon0, rad0;
    double maglat0, maglon0, magrad0;

    double dd;

    // Energy and velocity grids at output
    double v_tot_arr[NUM_E], E_tot_arr[NUM_E];


    rayT* r_cur;    // Current interpolated ray (single timestep)
    double* start_pos;

    rayF* ray;
    
    // Array of pointers to current rays of interest
    // (does this copy or just point? I hope it points.)
    rayF* cur_rays[8];

    // FILE * outputFile;
    

    // Default parameters:
    inpFileName = "input.ray";
    outFileName = "output.ray";
    dumpFileName= "fieldline_dump_igrf_tsyg.dat";

    int iyr;     // Year
    int idoy;      // Day of year
    double isec; // Seconds into day (UTC)

    // Fine-scale interpolation steps:
    int num_lats_fine  = 4;
    int num_lons_fine  = 4; 
    int num_freqs_fine = 4;

    double freq_step_size;
    double lat_step_size;
    double lon_step_size;
    
    EA_segment EA_array[NUM_EA];

    // Location to determine output at (geomagnetic)
    // To do: set this up as a grid
    double outLat = 60;
    double outLon = 180;

    int dump_field = 0;

     // Parse input arguments:
    int opt = 0;
    while ((opt = getopt(argc, argv, "i:o:t:u:v:a:b:c:d:")) != -1) {
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
            case 'd':
                dump_field = atoi(optarg);
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

    // Set up output grids (energy, velocity, space)

    //initialize the velocity and energy arrays
    for(int i=0; i<NUM_E; i++) {
        E_tot_arr[i] = pow(10, (E_EXP_BOT+ DE_EXP*i) ); // energy in eV
        v_tot_arr[i] = C*sqrt(1 - pow( (E_EL/(E_EL+E_tot_arr[i])) ,2) );
    }

    init_EA_array(EA_array, outLat, outLon, iyr, idoy, isec);


    int yearday = iyr*1000 + idoy;
    itime_in[0] = yearday;
    itime_in[1] = isec*1e3;


    // Dump the field line model:
    if (dump_field == 1) {
        int use_IGRF = 1;
        int use_tsyg = 0;
        int n_lats = 10;
        int n_lons = 4;


        dump_fieldlines(itime_in, n_lats, n_lons, use_IGRF, use_tsyg, dumpFileName);
    }





    // cout << "flash geo: ";
    // print_vector(vector<double>(flash_pos, flash_pos + 3));
    // double flash_I0 = 100e3;
    
    // int yearday = iyr*1000 + idoy;
    // itime_in[0] = yearday;
    // itime_in[1] = isec*1e3;

    // lat0 = D2R*flash_pos[1];
    // lon0 = D2R*flash_pos[2];
    // rad0 = flash_pos[0];

    // // Get flash position in SM coords:
    // pol_to_cart_d_(&lat0, &lon0, &rad0, tmp_coords);
    // mag_to_sm_d_(itime_in, tmp_coords, flash_pos_sm);

    // cout << "SM (libxformd): ";
    // print_vector(vector<double>(flash_pos_sm, flash_pos_sm + 3));








    // // Load the rayfile:
    // raylist = read_rayfile(inpFileName);


    // // Preprocess ray files: 
    // for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter){
    
    //     // cout << "Ray origin (SM):\n";

    //     ray = &(iter->second);

    //     start_pos = &(ray->pos[0].data()[0]);
    //     // print_vector(ray->pos[0]);


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

    //     // cout << "dd: " << ray->inp_pwr << "\n";

    //     // Calculate Stix parameters:
    //     calc_stix_parameters(ray);
    // }



    // // Choose the 8 corner rays we'll work with (this is where you'll iterate
    // // over the larger set)
    // cur_rays[0] = &(raylist.at(1));
    // cur_rays[1] = &(raylist.at(2));
    // cur_rays[2] = &(raylist.at(3));
    // cur_rays[3] = &(raylist.at(4));
    // cur_rays[4] = &(raylist.at(5));
    // cur_rays[5] = &(raylist.at(6));
    // cur_rays[6] = &(raylist.at(7));
    // cur_rays[7] = &(raylist.at(8));


    // cout << "sanity check: " << cur_rays[1]->stixB[0] << "\n";


    // // for (int dd = 0; dd < 8; dd++) {
    // //     cout << "Stix Params [0]: ";
    // //     cout << cur_rays[dd]->stixR[1] << " ";
    // //     cout << cur_rays[dd]->stixL[1] << " ";
    // //     cout << cur_rays[dd]->stixP[1] << " ";
    // //     cout << cur_rays[dd]->stixS[1] << " ";
    // //     cout << cur_rays[dd]->stixD[1] << " ";
    // //     cout << cur_rays[dd]->stixA[1] << " ";
    // //     cout << cur_rays[dd]->stixB[1] << "\n";
    // //     // cout << cur_rays[dd]->time.size() << " ";
    // // }

    // cout << "\n";

    // // Get the length of the shortest ray in the batch:
    // int tmaxes[] = {cur_rays[0]->time.size(),
    //                 cur_rays[1]->time.size(),
    //                 cur_rays[2]->time.size(),
    //                 cur_rays[3]->time.size(),
    //                 cur_rays[4]->time.size(),
    //                 cur_rays[5]->time.size(),
    //                 cur_rays[6]->time.size(),
    //                 cur_rays[7]->time.size()};
    // int tmax = *min_element(tmaxes, tmaxes + 8);

    // cout << "tmax is: " << tmax << "\n";

    // // Determine fine-scale grid steps:
    //     // (do this later)


    // for (int tt = 0; tt < 1; tt++) {
    //     cout << "t = " << tt << "\n";
    //     // Interpolate on fine-scale grid:
    //     for (double ii=0; ii <= 1; ii+=1./num_freqs_fine) {         // freqs
    //         for (double jj=0; jj <= 1; jj+= 1./num_lats_fine) {     // lats
    //             for (double kk=0; kk <= 1; kk+= 1./num_lons_fine) { // lons

    //                 // pos_interp = {0.,0.,0.};
                    
    //                 r_cur = new rayT;
    //                 // weight_ind = {ii, jj, kk};
    //                 int t_ind = 0;

    //                 interp_ray_fine(cur_rays, ii, jj, kk, t_ind, r_cur);

    //                 // Get magnetic lat:
    //                 sm_to_mag_d_(itime_in, r_cur->pos, tmp_coords2);
    //                 cart_to_pol_d_(tmp_coords2, &maglat0, &maglon0, &magrad0);
    //                 maglat0 = R2D*maglat0; maglon0 = R2D*maglon0;

    //                 // Print for debugging:
    //                 // printf("(%g, %g, %g)\tmag lat, lon: (%g, %g): freq: %g hz\n",ii, jj, kk, maglat0,maglon0, (r_cur->w)/(2.*PI));

    //             }   // kk
    //         }   // jj
    //     }   // ii
    // } // tt












    // for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter){

    //     cout << "debuggingd: " << iter->second.inp_pwr << "\n";
    // }
    return 0; // Return statement.
} // Closing Main.