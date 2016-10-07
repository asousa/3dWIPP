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

        // iyr = 2010; idoy = 1; isec = 3600;


    // cout << "---- Input Parameters ----\n";
    // cout << "input file: " << inpFileName << "\n";
    // cout << "output file: " << outFileName << "\n";
    // cout << "year: " << iyr << "\n";
    // cout << "day: " << idoy << "\n";
    // cout << "sec of day: " << isec << "\n";
    // cout << "---- 3D WIPP ----\n";

    // tmp_in[0] =-3653337.8925629;
    // tmp_in[1] =-2740585.99220064;
    // tmp_in[2] =4452100.02831447;
    
          
    // cout << "x_in: " << x_in[0] << ", " << x_in[1] << ", " << x_in[2] << "\n";

    // sm_to_geo_d_(itime_in, x_in, tmp_coords);
    // cart_to_pol_d_(tmp_coords, &lat0, &lon0, &rad0);
    // printf("ray origin (xformd): %g, %g, %g\n",lat0*R2D, lon0*R2D, rad0/1000.);

    // int iyr = int(floor(itime_in[0]/1000.0));
    // int idoy= int(fmod(itime_in[0], float(iyr)));
    // double secs = itime_in[1]*1e-3;

    // cout << "iyr: " << iyr << ", idoy: " << idoy << "\n";
    // sm2geo1_(&iyr, &idoy, &isec, x_in, x_out);
    // car_sph_(x_out, &rad0, &lat0, &lon0);
    // // cart_to_pol_d_(x_out, &lat0, &lon0, &rad0);
    // printf("ray origin (other thing): %g, %g, %g\n",lat0, lon0, rad0/1000.);

    // double in_lat = flash_pos[1];
    // double in_lon = flash_pos[2];
    // double in_rad = flash_pos[0];

    // double in_coords = {in_rad, in_lat, in_lon};


    // flash_pos = {1, 45, 17};


    cout << "flash geo: ";
    print_vector(vector<double>(flash_pos, flash_pos + 3));
    double flash_I0 = 100e3;
    
    ntimes=1; c1 = 8; c2 = 4;
    coord_trans_vec1_(&ntimes, &c1, &c2, &iyr, &idoy, &isec, flash_pos, flash_pos_sm);
    

    // cout << "Geo cartesian: " << tmp_coords[0] << " " << tmp_coords[1] << " " << tmp_coords[2] << "\n";
    // geo2sm1_(&iyr, &idoy, &isec, tmp_coords, x);
    cout << "flash SM: ";
    print_vector(vector<double>(flash_pos_sm, flash_pos_sm + 3));

    ntimes=1; c1 = 4; c2 = 8;
    coord_trans_vec1_(&ntimes, &c1, &c2, &iyr, &idoy, &isec, flash_pos_sm, flash_pos);

    cout << "flash geo (back): ";
    print_vector(vector<double>(flash_pos, flash_pos + 3));



    int yearday = iyr*1000 + idoy;
    // cout << yearday << "\n";

    itime_in[0] = yearday;
    itime_in[1] = isec;

    lat0 = D2R*flash_pos[1];
    lon0 = D2R*flash_pos[2];
    rad0 = flash_pos[0];

    pol_to_cart_d_(&lat0, &lon0, &rad0, tmp_coords);
    geo_to_sm_d_(itime_in, tmp_coords, tmp_coords2);

    // cout << "geo cart: "; 
    // print_vector(vector<double>(tmp_coords, tmp_coords + 3));

    cout << "SM (libxformd): ";
    print_vector(vector<double>(tmp_coords2, tmp_coords2 + 3));



    // cout << "SM cartesian: " << tmp_coords2[0] << " " << flash_pos_sm[1] << " " << flash_pos_sm[2] << "\n";



    // cout << " ----- Coordinate transform idiot check: \n";
    // long ntimes, c1, c2;

    // // Geo spherical inputs: Re, lat, lon:

    // x_in = {1, 1, 0};
    // cout << "X GEO Sph: " << x_in[0] << " " << x_in[1] << " " << x_in[2] << "\n";

    // // Coordinate transformation:
    // ntimes=1; c1 = 7; c2 = 4;
    // cout << "Conversion mode: " << c1 << " -> " << c2 << "\n";
    // coord_trans_vec1_(&ntimes, &c1, &c2, &iyr, &idoy, &isec, x_in, x_out);

    // cout << "X SM Cart: " << x_out[0] << " " << x_out[1] << " " << x_out[2] << "\n";
    
    // double y_in[3];
    // double y_out[3];
    // y_in = {1, 0, 0};
    // cout << "Y GEO Sph: " << y_in[0] << " " << y_in[1] << " " << y_in[2] << "\n";

    // // Coordinate transformation:
    // // ntimes=1; c1 = 7; c2 = 6;
    // coord_trans_vec1_(&ntimes, &c1, &c2, &iyr, &idoy, &isec, y_in, y_out);

    // cout << "Y SM Cart: " << y_out[0] << " " << y_out[1] << " " << y_out[2] << "\n";

    // // Find angle between the two vectors, in SM coordinates:

    // Vector3d v1 = Map<VectorXd>(x_out,3,1);
    // Vector3d v2 = Map<VectorXd>(y_out,3,1);
    
    // double theta = R2D*acos(v1.dot(v2)/(v1.norm()*v2.norm()));
    // cout << "theta: " << theta << " deg\n";



    // sph_car_(&in_rad, &in_lat, &in_lon, tmp_coords);
    // cout << "geo cartesian: " << tmp_coords[0] << " " << tmp_coords[1] << " " << tmp_coords[2] << "\n";





    // // Load the rayfile:
    // raylist = read_rayfile(inpFileName);

    // vector <double> start_pos;
    // // double start_pos_geo[3];

    // // Parse rays, make a nice list of the start coords and frequencies
    // for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter){
        
    //     // start_pos = iter->second.pos[0].front();

    //     ntimes=1; c1 = 4; c2 = 7;
    //     cout << "Conversion mode: " << c1 << " -> " << c2 << "\n";
    //     coord_trans_vec1_(&ntimes, &c1, &c2, &iyr, &idoy, &isec, iter->second.pos[0].data(), tmp_coords);
    //     ntimes=1; c1 = 4; c2 = 6;
    //     coord_trans_vec1_(&ntimes, &c1, &c2, &iyr, &idoy, &isec, iter->second.pos[0].data(), tmp_coords2);
    //     car_sph_(tmp_coords2, &magrad0, &maglat0, &maglon0);

    //     printf("SM ray: %g, %g, %g\n",iter->second.pos[0][0],iter->second.pos[0][1],iter->second.pos[0][2]);

    //     // sm2geo1_(&iyr, &idoy, &isec, &(iter->second.pos[0][0]), tmp_coords);
    //     // car_sph_(tmp_coords, &rad0, &lat0, &lon0);
    //     // geo2mag1_(&iyr, tmp_coords, tmp_coords2);
    //     // car_sph_(tmp_coords2, &magrad0, &maglat0, &maglon0);

    //     printf("GEO lat: %g lon: %g alt: %g\n",tmp_coords[1], tmp_coords[2], tmp_coords[0]);
    //     printf("MAG lat: %g lon: %g alt: %g\n",maglat0,maglon0,magrad0);
        
    //     Vector3d v1 = Map<VectorXd>(iter->second.pos[0].data(),3,1);
    //     Vector3d v2 = Map<VectorXd>(flash_pos_sm,3,1);
        
    //     double theta = R2D*acos(v1.dot(v2)/(v1.norm()*v2.norm()));
    //     cout << "theta: " << theta << " deg\n";


    //     // dd = input_power_scaling(flash_pos_sm, start_pos, maglat0, iter->second.w, flash_I0);

    // }



    return 0; // Return statement.
} // Closing Main.