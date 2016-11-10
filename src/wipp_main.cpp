#include <wipp.h>

using namespace std;
using namespace Eigen;


int main(int argc, char *argv[]) 
{
    map <int, rayF> raylist_hi;
    map <int, rayF> raylist_low;

    double flash_pos[3];  // Flash position, geographic spherical coords (r, lat, lon)
    double flash_pos_sm[3];
    double x_in[3];
    double x_out[3];

    string ray_inp_dir;
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

    map <pair<int, int>, cellT> crossing_db[NUM_EA];
    pair<int, int> grid_ind;
    int t_grid, f_grid;


    map <int, vector<double> > start_locs;
    vector < vector<int> > adjacency_list;  // ( n x 4 ) list of adjacent ray indexes

    // Array of pointers to current rays of interest
    // (does this copy or just point? I hope it points.)
    rayF* cur_rays[8];

    // Array of pointers to single-timestep frames
    rayT cur_frames[8];
    rayT prev_frames[8];

    // FILE * outputFile;
    

    // Default parameters:
    ray_inp_dir = "/shared/users/asousa/WIPP/3dWIPP/outputs/";
    // inpFileName = "input.ray";
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

    // double freq_step_size =10;    // Hz
    // double lat_step_size = 10;    // km
    // double lon_step_size = 100;   // km
    
    EA_segment EA_array[NUM_EA];

    // Location to determine output at (geomagnetic)
    // To do: set this up as a grid
    double outLat = 50;
    double outLon = 2;
    int model_number = 0; // Magnetic field model

    int dump_field = 0;

     // Parse input arguments:
    int opt = 0;
    while ((opt = getopt(argc, argv, "i:o:t:u:v:a:b:c:d:")) != -1) {
        switch(opt) {
            case 'i':
            // input filename:
            // cout << "inp filename!\n";
                ray_inp_dir = (string) optarg;
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
    cout << "ray input directory: " << ray_inp_dir << "\n";
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
    double flash_I0 = -100e3;
    
    lat0 = D2R*flash_pos[1];
    lon0 = D2R*flash_pos[2];
    rad0 = flash_pos[0];

    // Get flash position in SM coords:
    pol_to_cart_d_(&lat0, &lon0, &rad0, tmp_coords);
    mag_to_sm_d_(itime_in, tmp_coords, flash_pos_sm);

    cout << "SM (libxformd): ";
    print_vector(vector<double>(flash_pos_sm, flash_pos_sm + 3));


    double input_freqs[] = FREQ_VEC;  // Currently set in const.h


    cout << "input frequencies: "; print_array(input_freqs, 2);




    // Load the rayfile:
    // raylist_low = read_rayfile(inpFileName.str());
    // raylist_hi  = read_rayfile(inpFileName.str());

    // for (int f_ind = 0; f_ind < NUM_FREQS; f_ind++) {

    // Load upper frequency rays ------------------------------------------:
    ostringstream inpFileName;
    inpFileName << ray_inp_dir << "/rayout_" << input_freqs[1] << "_damped.ray";
    raylist_hi  = read_rayfile(inpFileName.str());
    
    // inpFileName << ray_inp_dir << "/rayout_" << input_freqs[f_ind - 1] << "_damped.ray";
    // raylist_low  = read_rayfile(inpFileName.str());

    // Preprocess ray files: 
    for(map<int,rayF>::iterator iter = raylist_hi.begin(); iter != raylist_hi.end(); ++iter){
        ray = &(iter->second);
        start_pos = &(ray->pos[0].data()[0]);

        // Get starting coordinates in geomagnetic:
        sm_to_mag_d_(itime_in, start_pos, tmp_coords2);
        cardeg(tmp_coords2);

        // printf("MAG lat: %g lon: %g alt: %g\n",tmp_coords2[1],tmp_coords2[2],tmp_coords2[0]);
        start_locs.insert(make_pair(iter->first, vector<double>(tmp_coords2, tmp_coords2 + 3)));
        iter->second.in_radius = tmp_coords2[0]; iter->second.in_lat = tmp_coords2[1]; iter->second.in_lon = tmp_coords2[2];

        // Calculate Stix parameters:
        calc_stix_parameters(ray);
    }

    // Load lower frequency rays -------------------------------------------:
    inpFileName.str(""); inpFileName.clear();

    inpFileName << ray_inp_dir << "/rayout_" << input_freqs[0] << "_damped.ray";
    raylist_low  = read_rayfile(inpFileName.str());
    
    // inpFileName << ray_inp_dir << "/rayout_" << input_freqs[f_ind - 1] << "_damped.ray";
    // raylist_low  = read_rayfile(inpFileName.str());

    // Preprocess ray files: 
    for(map<int,rayF>::iterator iter = raylist_low.begin(); iter != raylist_low.end(); ++iter){
        ray = &(iter->second);
        start_pos = &(ray->pos[0].data()[0]);

        // Get starting coordinates in geomagnetic:
        sm_to_mag_d_(itime_in, start_pos, tmp_coords2);
        cardeg(tmp_coords2);

        printf("MAG lat: %g lon: %g alt: %g\n",tmp_coords2[1],tmp_coords2[2],tmp_coords2[0]);
        ray->in_radius = tmp_coords2[0]; ray->in_lat = tmp_coords2[1]; ray->in_lon = tmp_coords2[2];
        start_locs.insert(make_pair(iter->first, vector<double>(tmp_coords2, tmp_coords2 + 3)));

        printf("(from ray): %g %g %g\n",ray->in_radius, ray->in_lat, ray->in_lon);
        // Calculate Stix parameters:
        calc_stix_parameters(ray);
    }


    // Find all sets of adjacent rays to iterate over:
    adjacency_list = find_adjacent_rays(start_locs);

    cout << "Found " << adjacency_list.size() << " sets of adjacent guide rays\n";

    // Iterate over each set of adjacent guide rays:
    for (int adj_row =0; adj_row < adjacency_list.size(); adj_row++) {
        // Choose the 8 corner rays we'll work with (this is where you'll iterate
        // over the larger set)
        cur_rays[0] = &(raylist_low.at(adjacency_list[adj_row][0]));
        cur_rays[1] = &(raylist_low.at(adjacency_list[adj_row][1]));
        cur_rays[2] = &(raylist_low.at(adjacency_list[adj_row][2]));
        cur_rays[3] = &(raylist_low.at(adjacency_list[adj_row][3]));
        cur_rays[4] = &(raylist_hi.at( adjacency_list[adj_row][0]));
        cur_rays[5] = &(raylist_hi.at( adjacency_list[adj_row][1]));
        cur_rays[6] = &(raylist_hi.at( adjacency_list[adj_row][2]));
        cur_rays[7] = &(raylist_hi.at( adjacency_list[adj_row][3]));
        
        // Find minimum and maximum frequencies, start lats, and start lons:
        wmin   = cur_rays[0]->w;         wmax = cur_rays[0]->w;
        latmin = cur_rays[0]->in_lat;  latmax = cur_rays[0]->in_lat;
        lonmin = cur_rays[0]->in_lon;  lonmax = cur_rays[0]->in_lon;
        tmax   = cur_rays[0]->time.back();
        double in_lat, in_lon;
        for (int i=1; i < 8; i++) {
            in_lon = cur_rays[i]->in_lon;
            in_lat = cur_rays[i]->in_lat;
            // cout << "HEY FUCKER " << in_lon << " " << in_lat << "\n"; 
            if (in_lon >= 360) { in_lon -= 360; }

            if (cur_rays[i]->w < wmin)  { wmin = cur_rays[i]->w; } 
            if (cur_rays[i]->w > wmax ) { wmax = cur_rays[i]->w; }
            if (cur_rays[i]->in_lat < latmin ) { latmin = cur_rays[i]->in_lat; }
            if (cur_rays[i]->in_lat > latmax ) { latmax = cur_rays[i]->in_lat; }
            if (in_lon < lonmin ) { lonmin = in_lon; }
            if (in_lon > lonmax ) { lonmax = in_lon; }
            if (cur_rays[i]->time.back() < tmax) { tmax = cur_rays[i]->time.back(); }

        }

        // starting separation in lat, lon directions (meters)
        dlat = D2R*R_E*(latmax - latmin);
        dlon = D2R*R_E*(lonmax - lonmin)*cos(D2R*(latmax + latmin)/2.);
        dw   = wmax - wmin;

        cout << "lon: " << lonmax << ", " << lonmin << "\n";
        cout << "lat: " << latmax << ", " << latmin << "\n";
        cout << "dlat: " << dlat << " dlon: " << dlon << "\n";
        cout << "dw: " << dw << "\n";

        // Scale the input power by dlat, dlon, dw:
        // (ray spacing may not be consistent)
        for (int i=1; i < 8; i++) {
            cur_rays[i]->inp_pwr = input_power_scaling(flash_pos_sm, cur_rays[i]->pos[0].data(),
                                   cur_rays[i]->in_lat, cur_rays[i]->w, flash_I0);

            cur_rays[i]->inp_pwr *= (dlat)*(dw/(2*PI)*0.877);       // This matches Jacob's power scaling (but I disagree with it)
            cout << "in pwr (post scale): " << cur_rays[i]->inp_pwr << "\n";
        }


        // // Always do at least 2 steps in each axis (corner rays)
        num_freqs_fine = max(2, (int)floor( (wmax - wmin)/(2*PI*FREQ_STEP_SIZE )));
        num_lats_fine  = max(2, (int)floor( (dlat*1e-3)/(LAT_STEP_SIZE) ));
        num_lons_fine  = max(1, (int)floor( (dlon*1e-3)/(LON_STEP_SIZE) ));


        // num_freqs_fine = 4;
        // num_lats_fine  = 10;
        // num_lons_fine  = 4;

        cout << "num steps: " << num_freqs_fine << ", " << num_lats_fine << ", " << num_lons_fine << "\n";

        crossing_log = fopen(crossingFileName.c_str(), "w");


        cout << "T_STEP: " << TIME_STEP << "\n";
        double hit_counter = 0;
        double crossing_counter = 0;

        cout << "checking for crossings...\n";

        
        cout << "input frequencies: ";

        for (int eh=0; eh < 8; eh++) {
            cout << cur_rays[eh]->w/(2*PI) << " ";
        }




        // Interpolate the first frames:
        for (int zz=0; zz<8; zz++) { interp_rayF(cur_rays[zz], &(prev_frames[zz]), 0); }

        for (double tt=TIME_STEP; tt < tmax; tt+=TIME_STEP) {


        // for (int tt = 1; tt < tmax; tt++) {
            // interpolate current frames:
            for (int zz=0; zz<8; zz++) { interp_rayF(cur_rays[zz], &(cur_frames[zz]), tt); }

            // cout << "time: " << cur_frames[0].time << ", " << prev_frames[0].time << "\n";


            // Check each EA segment:
            for (int rr = 0; rr < NUM_EA; rr++) {
                // cout << "EA: " << EA_array[rr].lat << "\n";

                // Ignore anything that looks way out of range
                if (coarse_mask(cur_frames, prev_frames, EA_array[rr])) {
                // if(true){
                    // cout << "hit\n";
                    hit_counter ++;

                    // Interpolate on fine-scale grid:
                    for (double ii=0; ii <= 1; ii+=1./num_lons_fine) {         
                        for (double jj=0; jj < 1; jj+= 1./num_lats_fine) {     
                            for (double kk=0; kk < 1; kk+= 1./num_freqs_fine) { 
                                // cout << ii << " ";

                                // // Clear previous values
                                r_cur =  {};
                                r_prev = {};

                                // // (to do: Save r_curs to avoid having to recalculate it)
                                // interp_ray_positions(cur_rays, ii, jj, kk, tt,  &r_cur);
                                interp_ray_positions(cur_frames, ii, jj, kk, &r_cur);
                                interp_ray_positions(prev_frames,ii, jj, kk, &r_prev);



                                // // Bam -- we finally have some little rays to check for crossings.
                                if (crosses_EA(r_cur.pos, r_prev.pos, EA_array[rr])) {
                                    crossing_counter++;

                                    interp_ray_data(cur_frames, ii, jj, kk, &r_cur);
                                    interp_ray_data(prev_frames,ii, jj, kk, &r_prev);

                                    
                                    
                                    // cout << "f_int: " << r_cur.w/(2*PI) << "\n";
                                    fprintf(crossing_log, "%g %g %g %g %g %g\n",
                                        r_cur.pos[0], r_cur.pos[1], r_cur.pos[2], 
                                        r_prev.pos[0], r_prev.pos[1], r_prev.pos[2]);

                                    // store time and frequency for the middle of this interpolation
                                    r_cur.dt = (r_cur.time - r_prev.time);
                                    r_cur.dlat = dlat;
                                    r_cur.dlon = dlon;
                                    r_cur.ds   = (r_cur.pos - r_prev.pos).norm()*R_E;   // Jacob uses ds between the EA segments... hm
                                    
                                    // cout << "tt: " << tt << "t_cur: " << r_cur.time << "\n";
                                    t_grid = floor(r_cur.time/(TIME_STEP));
                                    f_grid = floor(kk*num_freqs_fine); //floor(r_cur.w/num_freqs_fine);
                                                                        // cout << "t_grid: " << t_grid << " f_grid: " << f_grid << "\n";
                                    grid_ind = make_pair(t_grid, f_grid);

                                    cellT cell_cur = new_cell(r_cur);

                                    cell_cur.Lsh = EA_array[rr].Lsh;
                                    cell_cur.lat = EA_array[rr].lat;
                                    if (crossing_db[rr].count(grid_ind)==0) {
                                        // If we haven't hit this same (time, freq, EA) combo yet, add it:
                                        crossing_db[rr].insert(make_pair(grid_ind, cell_cur));
                                    } else {
                                        // Else, sum the current frame with previous frames, so we can average:
                                        // add_rayT(&(crossing_db[rr].at(grid_ind)), &r_cur);
                                        add_cell(&(crossing_db[rr].at(grid_ind)), &cell_cur);
                                    }

                                }   // Crossings
                            }   // kk
                        }   // jj
                    }   // ii
                }   // Coarse mask
            }   // EA array

            // Step forward one frame:
            for (int zz=0; zz<8; zz++) { prev_frames[zz] = cur_frames[zz]; }
        } // tt




        cout << "hit counter: " << hit_counter << "\n";
        cout << "crossing counter: " << crossing_counter << "\n";

    } // Cur_rays

    


    // Calculate scattering at crossings:
    // cout << "Calculating resonances\n";
    for (int rr=NUM_EA-1; rr>=0; rr--) {

        printf("---------------- EA at lat: %2.2f ------------------\n", EA_array[rr].lat);
        for(map<pair<int,int>,cellT>::iterator iter = crossing_db[rr].begin(); iter != crossing_db[rr].end(); ++iter){
            // calc_resonance(&(iter->second), &(EA_array[rr]), v_tot_arr, da_N, da_S);
            calcRes(iter->second, da_N, da_S);
        }
    }




    // // Let's try this with the old crossings:
    // vector<cellT> listy = load_crossings(itime_in, "/shared/users/asousa/WIPP/WIPPv4/crossing_log.txt");
    // for (vector<cellT>::iterator cell = listy.begin(); cell!=listy.end(); cell++) {
    //     calcRes(*cell, da_N, da_S);
    // }


    // for (int i=0; i < NUM_TIMES; i++) {
    //     for (int j=0; j < NUM_E; j++) {
    //         if (da_N[j][i] > 0) {
    //             cout << i << " " << j << " " << sqrt(da_N[j][i]) << "\n";
    //         }  
    //     }  
    // }


    // Write pN, pS files:
    cout << "Saving pN, pS files\n";
    write_p_array(da_N, "pN.dat");
    write_p_array(da_S, "pS.dat");




    return 0; // Return statement.
} // Closing Main.