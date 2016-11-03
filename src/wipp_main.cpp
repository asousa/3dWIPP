#include <wipp.h>

using namespace std;
using namespace Eigen;


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

    // debugging and dump filenames
    string dumpFileName;
    string eaFileName;
    string crossingFileName;

    int itime_in[2];
    
    string itime_str;

    double tmp_coords[3];
    double tmp_coords2[3];

    double x[3];
    double lat0, lon0, rad0;
    double maglat0, maglon0, magrad0;

    double dd;

    double wmin, wmax, latmin, latmax, lonmin, lonmax;
    double dw, dlon, dlat;  
    int tmax;

    // Energy and velocity grids at output
    double v_tot_arr[NUM_E], E_tot_arr[NUM_E];


    rayT r_cur;    // Current interpolated ray (single timestep)
    rayT r_prev;   // Previous interpolated ray

    Vector3d l0, l1;
    double* start_pos;

    rayF* ray;
    
    // Output arrays (change in pitch angle, north and south hemispheres)
    double da_N[NUM_E][NUM_TIMES] = {0};
    double da_S[NUM_E][NUM_TIMES] = {0};


    // Array of pointers to current rays of interest
    // (does this copy or just point? I hope it points.)
    rayF* cur_rays[8];

    // FILE * outputFile;
    

    // Default parameters:
    inpFileName = "input.ray";
    outFileName = "output.ray";
    dumpFileName= "fieldline_dump.dat";
    eaFileName  = "ea_dump.dat";
    crossingFileName = "crossing_log.txt";

    FILE* crossing_log;

    int iyr;     // Year
    int idoy;      // Day of year
    double isec; // Seconds into day (UTC)

    // Fine-scale interpolation steps:
    int num_lats_fine;
    int num_lons_fine;   
    int num_freqs_fine; 

    double freq_step_size = 10;   // Hz
    double lat_step_size = 100;   // km
    double lon_step_size = 100;   // km
    
    EA_segment EA_array[NUM_EA];

    // Location to determine output at (geomagnetic)
    // To do: set this up as a grid
    double outLat = 50;
    double outLon = 0;
    int model_number = 0; // Magnetic field model

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


    int yearday = iyr*1000 + idoy;
    itime_in[0] = yearday;
    itime_in[1] = isec*1e3;


    // Set up output grids (energy, velocity, space)

    //initialize the velocity and energy arrays
    for(int i=0; i<NUM_E; i++) {
        E_tot_arr[i] = pow(10, (E_EXP_BOT+ DE_EXP*i) ); // energy in eV
        v_tot_arr[i] = C*sqrt(1 - pow( (E_EL/(E_EL+E_tot_arr[i])) ,2) );
    }

    init_EA_array(EA_array, outLat, outLon, itime_in, model_number);
    dump_EA_array(EA_array, eaFileName);

    // Dump the field line model:
    if (dump_field == 1) {
        int n_lats = 10;
        int n_lons = 4;

        dump_fieldlines(itime_in, n_lats, n_lons, model_number, dumpFileName);
    }

    srand(time(0));



    cout << "flash geo: ";
    print_vector(vector<double>(flash_pos, flash_pos + 3));
    double flash_I0 = 100e3;
    
    // int yearday = iyr*1000 + idoy;
    // itime_in[0] = yearday;
    // itime_in[1] = isec*1e3;

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


    // Preprocess ray files: 
    for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter){
    
        ray = &(iter->second);

        start_pos = &(ray->pos[0].data()[0]);

        // Get magnetic lat:
        sm_to_mag_d_(itime_in, start_pos, tmp_coords2);
        cart_to_pol_d_(tmp_coords2, &maglat0, &maglon0, &magrad0);
        maglat0 = R2D*maglat0; maglon0 = R2D*maglon0;
        printf("MAG lat: %g lon: %g alt: %g\n",maglat0,maglon0,magrad0);
        
        // Get power scaling:
        // (still need to multiply by space + freq bin sizes)
        // ray->inp_pwr = input_power_scaling(flash_pos_sm, start_pos, maglat0, iter->second.w, flash_I0);
        ray->in_radius = magrad0;
        ray->in_lat = maglat0;
        ray->in_lon = maglon0;

        // Calculate Stix parameters:
        calc_stix_parameters(ray);
    }


    // Choose the 8 corner rays we'll work with (this is where you'll iterate
    // over the larger set)
    cur_rays[0] = &(raylist.at(1));
    cur_rays[1] = &(raylist.at(2));
    cur_rays[2] = &(raylist.at(3));
    cur_rays[3] = &(raylist.at(4));
    cur_rays[4] = &(raylist.at(5));
    cur_rays[5] = &(raylist.at(6));
    cur_rays[6] = &(raylist.at(7));
    cur_rays[7] = &(raylist.at(8));

    cout << "\n";



    double tsyg_params[10] = {0};
    double VG[3];
   
    // int model_number = 0;

    //  // Setup for IGRF:    
    // load_TS05_params(itime_in, tsyg_params, VG);
    // init_igrf(itime_in);

    // // set DST:
    // tsyg_params[1] = -20;



    // double btmp[3];
    // // Confirm we're using the same B model:
    // for (int i =0; i < 8; i++) {
    //     for (int t=0; t < cur_rays[i]->time.size(); t++) {

    //         bmodel(itime_in, cur_rays[i]->pos[t].data(), tsyg_params, model_number, btmp);

    //         Vector3d nt = Map<Vector3d>(cur_rays[i]->n[t].data(),3,1);
    //         Vector3d Bt = Map<Vector3d>(cur_rays[i]->B0[t].data(),3,1);
    //         double psi = -90*D2R - nt.dot(Bt)/(nt.norm()*Bt.norm());
    //         cout << "t= " << cur_rays[i]->time[t] << " psi: " << R2D*psi << "\n";
    //         // cout << "b(model) ";
    //         // print_array(btmp,3);
    //         // cout << "b(raytracer) ";
    //         // print_array(cur_rays[i]->B0[t].data(), 3);
    //     }
    // }





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


    // Find minimum and maximum frequencies, start lats, and start lons:
    wmin   = cur_rays[0]->w;         wmax = cur_rays[0]->w;
    latmin = cur_rays[0]->in_lat;  latmax = cur_rays[0]->in_lat;
    lonmin = cur_rays[0]->in_lon;  lonmax = cur_rays[0]->in_lon;
    tmax   = cur_rays[0]->time.size();
    double in_lat, in_lon;
    for (int i=1; i < 8; i++) {
        in_lon = cur_rays[i]->in_lon;
        if (in_lon >= 360) { in_lon -= 360; }

        if (cur_rays[i]->w < wmin)  { wmin = cur_rays[i]->w; } 
        if (cur_rays[i]->w > wmax ) { wmax = cur_rays[i]->w; }
        if (cur_rays[i]->in_lat < latmin ) { latmin = cur_rays[i]->in_lat; }
        if (cur_rays[i]->in_lat > latmax ) { latmax = cur_rays[i]->in_lat; }
        if (in_lon < lonmin ) { lonmin = in_lon; }
        if (in_lon > lonmax ) { lonmax = in_lon; }
        if (cur_rays[i]->time.size() < tmax) { tmax = cur_rays[i]->time.size(); }
    }

    // starting separation in lat, lon directions (meters)
    dlat = D2R*R_E*(latmax - latmin);
    dlon = D2R*R_E*(lonmax - lonmin)*cos(D2R*(latmax + latmin)/2.);
    dw   = wmax - wmin;

    cout << "lon: " << lonmax << ", " << lonmin << "\n";


    // Scale the input power by dlat, dlon, dw:
    // (ray spacing may not be consistent)
    for (int i=1; i < 8; i++) {
        cur_rays[i]->inp_pwr = input_power_scaling(flash_pos_sm, cur_rays[i]->pos[0].data(),
                               cur_rays[i]->in_lat, cur_rays[i]->w, flash_I0);
        // cout << "in pwr (pre scale): " << cur_rays[i]->inp_pwr << "\n";
        // cur_rays[i]->inp_pwr *= (dlat*dlon)*(2*PI*dw);
        cur_rays[i]->inp_pwr *= (2*PI*dw);

        // cout << "in pwr (post scale): " << cur_rays[i]->inp_pwr << "\n";
    }


    // Always do at least 2 steps in each axis (corner rays)
    num_freqs_fine = max(2, (int)ceil( (wmax - wmin)/(2*PI*freq_step_size )));
    num_lats_fine  = max(2, (int)ceil( (dlat*1e-3)/(lat_step_size) ));
    num_lons_fine  = max(2, (int)ceil( (dlon*1e-3)/(lon_step_size) ));

    cout << "num steps: " << num_freqs_fine << ", " << num_lats_fine << ", " << num_lons_fine << "\n";

    crossing_log = fopen(crossingFileName.c_str(), "w");



    double hit_counter = 0;
    double crossing_counter = 0;

    // Check each EA segment:
    for (int rr = 0; rr < NUM_EA; rr++) {
        cout << "EA: " << EA_array[rr].lat << "\n";

        for (int tt = 1; tt < tmax; tt++) {
            // cout << "t = " << tt << "\n";

            
            // Ignore anything that looks way out of range
            if (coarse_mask(cur_rays, tt, EA_array[rr])) {

                hit_counter ++;

                // Interpolate on fine-scale grid:
                for (double ii=0; ii <= 1; ii+=1./num_freqs_fine) {         // freqs
                    for (double jj=0; jj <= 1; jj+= 1./num_lats_fine) {     // lats
                        for (double kk=0; kk <= 1; kk+= 1./num_lons_fine) { // lons

                            // Clear previous values
                            r_cur = {};
                            r_prev = {};
                            // (to do: Save r_curs to avoid having to recalculate it)
                            interp_ray_fine(cur_rays, ii, jj, kk, tt,  &r_cur);
                            interp_ray_fine(cur_rays, ii, jj, kk, tt-1,&r_prev);

                            // l0 = r_cur.pos;// Map<VectorXd>(r_cur.pos, 3,  1);
                            // l1 = r_prev.pos;//Map<VectorXd>(r_prev.pos, 3, 1);

                            // Bam -- we finally have some little rays to check for crossings.
                            if (crosses_EA(r_cur.pos, r_prev.pos, EA_array[rr])) {
                                crossing_counter++;
                                // cout << l0.transpose() << " " << l1.transpose() << "\n";
                                fprintf(crossing_log, "%g %g %g %g %g %g\n",
                                    r_cur.pos[0], r_cur.pos[1], r_cur.pos[2], 
                                    r_prev.pos[0], r_prev.pos[1], r_prev.pos[2]);

                                // store time and frequency for the middle of this interpolation
                                r_cur.dt = (r_cur.time - r_prev.time);
                                r_cur.dw = (r_cur.w - r_prev.w);
                                r_cur.dlat = dlat;
                                r_cur.dlon = dlon;

                                calc_resonance(&r_cur, &(EA_array[rr]), v_tot_arr, da_N, da_S);

                            }   // Crossings
                        }   // kk
                    }   // jj
                }   // ii
            }   // Coarse mask
        }   // EA array
    } // tt


    cout << "hit counter: " << hit_counter << "\n";
    cout << "crossing counter: " << crossing_counter << "\n";


    // Write pN, pS files:

    write_p_array(da_N, "pN.dat");
    write_p_array(da_S, "pS.dat");



    // for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter){

    //     cout << "debuggingd: " << iter->second.inp_pwr << "\n";
    // }
    return 0; // Return statement.
} // Closing Main.