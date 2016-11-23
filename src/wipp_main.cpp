#include <wipp.h>

using namespace std;
using namespace Eigen;

#pragma omp 
int main(int argc, char *argv[]) 
{
    map <int, rayF> raylist_hi;
    map <int, rayF> raylist_low;

    double flash_pos[3];  // Flash position, geographic spherical coords (r, lat, lon)
    double flash_pos_sm[3];
    double x_in[3];
    double x_out[3];

    double flash_I0;

    string ray_inp_dir;
    string out_dir;

    // debugging and dump filenames
    string dumpFileName;
    ostringstream eaFileName;
    string crossingFileName;

    string low_file, high_file;

    int itime_in[2];
    
    double tmp_coords[3];
    double tmp_coords2[3];

    double x[3];
    double lat0, lon0, rad0;
    double maglat0, maglon0, magrad0;

    double dd;

    double wmin, wmax, latmin, latmax, lonmin, lonmax;
    double dw, dlon, dlat;  
    int tmax;

 
    rayT r_cur;    // Current interpolated ray (single timestep)
    rayT r_prev;   // Previous interpolated ray

    // Vector3d l0, l1;
    double* start_pos;

    rayF* ray;
    
    // // Output arrays (change in pitch angle, north and south hemispheres)
    double da_N[NUM_E][NUM_TIMES] = {0};
    double da_S[NUM_E][NUM_TIMES] = {0};

    map <pair<int, int>, cellT> crossing_db[NUM_EA];
    pair<int, int> grid_ind;
    int t_grid, f_grid;

    map <int, vector<double> > start_locs;
    vector < vector<int> > adjacency_list;  // ( n x 4 ) list of adjacent ray indexes

    // Array of pointers to current rays of interest
    rayF* cur_rays[8];

    // Array of pointers to single-timestep frames
    rayT cur_frames[8];
    rayT prev_frames[8];

    time_t run_tstart, run_tend;

    double frame_area, initial_area, cell_area, geometric_factor;

    int model_number = 0; // Magnetic field model
    int dump_field = 0;

    FILE* crossing_log;

    int iyr;       // Year
    int idoy;      // Day of year
    double isec;   // Seconds into day (UTC)

    // Fine-scale interpolation steps:
    int num_lats_fine;
    int num_lons_fine;   
    int num_freqs_fine; 
   
    vector<EA_segment> EA_array(NUM_EA);

    // Location to determine output at (geomagnetic)
    double out_lat = 50;
    double out_lon = 0;

    int opt = 0;
    int opt_index = 0;

    // Default parameters:
    out_dir = "/shared/users/asousa/WIPP/3dWIPP/outputs/";
    dumpFileName= "fieldline_dump.dat";
    crossingFileName = "crossing_log.txt";
    flash_pos = {1, 0, 0};
    flash_I0 = -100e3;
    out_lat = 45;
    out_lon = 0;
    iyr = 2010; idoy = 1; isec = 0;
    low_file  = "/shared/users/asousa/WIPP/3dWIPP/outputs/four_adjacent/rayout_1000_damped.ray";
    high_file = "/shared/users/asousa/WIPP/3dWIPP/outputs/four_adjacent/rayout_1100_damped.ray";



    // Parse input arguments:
    static struct option long_options[] =
    {
        {"out_dir",     required_argument,    0, 'a'},
        {"iyr",         required_argument,    0, 'b'},
        {"idoy",        required_argument,    0, 'c'},
        {"isec",        required_argument,    0, 'd'},
        {"f_alt",       required_argument,    0, 'e'},
        {"f_lat",       required_argument,    0, 'f'},
        {"f_lon",       required_argument,    0, 'g'},
        {"out_lat",     required_argument,    0, 'h'},
        {"out_lon",     required_argument,    0, 'i'},
        {"I0",          required_argument,    0, 'j'},
        {"low_file",    required_argument,    0, 'k'},
        {"hi_file" ,    required_argument,    0, 'l'},
        // These options set a flag:
        {"dump_field",  no_argument,    &dump_field, 1},
        {0, 0, 0, 0}
    };

    while (opt != -1) {
        opt = getopt_long (argc, argv, "a:b:c:d:e:f:g:h:i:j:k:", long_options, &opt_index);
        // cout << "opt is " << opt << "\n";
        switch(opt) {
            case 0:
            if (long_options[opt_index].flag != 0)      break;
            case 'a':   // out_dir
                out_dir = (string) optarg;              break;
            case 'b':   // iyr
                iyr = atoi(optarg);                     break;
            case 'c':   // idoy
                idoy= atoi(optarg);                     break;
            case 'd':   // isec
                isec= atoi(optarg);                     break;
            case 'e':   // flash altitude
                flash_pos[0] = strtod(optarg, NULL);    break;
            case 'f':   // flash latitude
                flash_pos[1] = strtod(optarg, NULL);    break;
            case 'g':   // flash longitude
                flash_pos[2] = strtod(optarg, NULL);    break;
            case 'h':   // out latitude
                out_lat = strtod(optarg, NULL);         break;
            case 'i':   // out longitude
                out_lon = strtod(optarg, NULL);         break;
            case 'j':
                flash_I0 = strtod(optarg, NULL);        break;
            case 'k':   // lower rayfile
                low_file= (string) optarg;              break;
            case 'l':   // upper rayfile
                high_file= (string) optarg;             break; 
            case '?':
                 printf("\nUnknown option: %s\n",opt);
            break;
        }
    }


    cout << "---- 3D WIPP ----\n\n";
    cout << "   Crossing Detection + pitch-angle deflection\n";
    cout << "---- Input Parameters ----\n";
    cout << "ray input directory: " << ray_inp_dir << "\n";
    cout << "output directory: " << out_dir << "\n";
    cout << "year: " << iyr << "\n";
    cout << "day: " << idoy << "\n";
    cout << "sec of day: " << isec << "\n";
    cout << "flash location: " << flash_pos[1] << ", " << flash_pos[2] << "\n";
    cout << "flash peak current: " << flash_I0*1e-3 << " kA\n";
    cout << "output location: " << out_lat << ", " << out_lon << "\n";
    cout << "\n";
    int yearday = iyr*1000 + idoy;
    itime_in[0] = yearday;
    itime_in[1] = isec*1e3;


// --------------- Set up output grid + EA array -----------------------
    EA_array = init_EA_array(out_lat, out_lon, itime_in, model_number);
    
    eaFileName.str(""); eaFileName.clear();
    eaFileName << out_dir << "/EA_dump_" << out_lat << "_" << out_lon << ".dat";
    dump_EA_array(EA_array, eaFileName.str());


    // Dump the field line model:
    if (dump_field == 1) {
        int n_lats = 10;
        int n_lons = 4;
        dump_fieldlines(itime_in, n_lats, n_lons, model_number, dumpFileName);
    }

// --------------- Get flash input coordinates -----------------------    
    lat0 = D2R*flash_pos[1];
    lon0 = D2R*flash_pos[2];
    rad0 = flash_pos[0];

    // Get flash position in SM coords:
    pol_to_cart_d_(&lat0, &lon0, &rad0, tmp_coords);
    mag_to_sm_d_(itime_in, tmp_coords, flash_pos_sm);


// Load upper frequency rays ------------------------------------------:
    raylist_hi  = read_rayfile(high_file);
    for(map<int,rayF>::iterator iter = raylist_hi.begin(); iter != raylist_hi.end(); ++iter){
        ray = &(iter->second);
        start_pos = &(ray->pos[0].data()[0]);

        // Get starting coordinates in geomagnetic:
        sm_to_mag_d_(itime_in, start_pos, tmp_coords2);
        cardeg(tmp_coords2);

        start_locs.insert(make_pair(iter->first, vector<double>(tmp_coords2, tmp_coords2 + 3)));
        ray->in_radius = tmp_coords2[0]; ray->in_lat = tmp_coords2[1]; ray->in_lon = tmp_coords2[2];

        // Calculate Stix parameters:
        calc_stix_parameters(ray);
    }

// Load lower frequency rays -------------------------------------------:
    raylist_low  = read_rayfile(low_file);
    for(map<int,rayF>::iterator iter = raylist_low.begin(); iter != raylist_low.end(); ++iter){
        ray = &(iter->second);
        start_pos = &(ray->pos[0].data()[0]);

        // Get starting coordinates in geomagnetic:
        sm_to_mag_d_(itime_in, start_pos, tmp_coords2);
        cardeg(tmp_coords2);

        ray->in_radius = tmp_coords2[0]; ray->in_lat = tmp_coords2[1]; ray->in_lon = tmp_coords2[2];

        // Calculate Stix parameters:
        calc_stix_parameters(ray);
    }

// Find all sets of adjacent rays to iterate over ----------------------:
    // (should be consistent between all frequency files... hmm... hmm...)
    adjacency_list = find_adjacent_rays(start_locs);
    cout << "Found " << adjacency_list.size() << " sets of adjacent guide rays\n";

    
    // Start with a fresh crossing db:
    for (int rr=0; rr < NUM_EA; ++rr) { crossing_db[rr].clear();}  

    // Iterate over each set of adjacent guide rays:
    for (int adj_row =0; adj_row < adjacency_list.size(); adj_row++) {
        // Choose the 8 corner rays to work with
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
        double avg_distance_from_flash = 0;

        for (int i=1; i < 8; i++) {
            in_lon = cur_rays[i]->in_lon;
            in_lat = cur_rays[i]->in_lat;

            avg_distance_from_flash += haversine_distance(in_lat, in_lon, flash_pos[1], flash_pos[2]);
            
            if (in_lon >= 360)                  { in_lon -= 360; }
            if (cur_rays[i]->w < wmin)          { wmin = cur_rays[i]->w; } 
            if (cur_rays[i]->w > wmax )         { wmax = cur_rays[i]->w; }
            if (cur_rays[i]->in_lat < latmin )  { latmin = cur_rays[i]->in_lat; }
            if (cur_rays[i]->in_lat > latmax )  { latmax = cur_rays[i]->in_lat; }
            if (in_lon < lonmin )               { lonmin = in_lon; }
            if (in_lon > lonmax )               { lonmax = in_lon; }
            if (cur_rays[i]->time.back() < tmax){ tmax = cur_rays[i]->time.back(); }
        }

        avg_distance_from_flash /= 8000.0;  // Average dist of the 8 corner rays, in km 


        // starting separation in lat, lon directions (meters)
        dlat = D2R*(R_E + H_IONO)*(latmax - latmin);
        dlon = D2R*(R_E + H_IONO)*(lonmax - lonmin)*cos(D2R*(latmax + latmin)/2.);
        dw   = wmax - wmin;

        // Only examine sets which are within our region of interest:
        if (avg_distance_from_flash <= MAX_GROUND_DISTANCE) {

            cout << "\n----- current rays: -----\n";
            cout << "lon: " << lonmax << ", " << lonmin << "\n";
            cout << "lat: " << latmax << ", " << latmin << "\n";
            cout << "f: " << wmax/(2*PI) << ", " << wmin/(2*PI) << "\n"; 
            cout << "Avg distance from flash: " << avg_distance_from_flash << "\n";


            // Scale the input power by dlat, dlon, dw:
            // (ray spacing may not be consistent)
            double inp_pwr = 0;

            for (int i=1; i < 8; i++) {
                cur_rays[i]->inp_pwr = input_power_scaling(flash_pos_sm, cur_rays[i]->pos[0].data(),
                                       cur_rays[i]->in_lat, cur_rays[i]->w, flash_I0);

                // cout << "inp_pwr: " << cur_rays[i]->inp_pwr << "\n";
                inp_pwr += cur_rays[i]->inp_pwr;

                // This matches Jacob's power scaling
                // cur_rays[i]->inp_pwr *= (dlat)*(dw/(2*PI)*0.877);      
                // cout << "in pwr (post scale): " << cur_rays[i]->inp_pwr << "\n";
            }

            // double tot_pwr = total_input_power(flash_pos_sm, flash_I0, 
            //                             latmin, latmax, lonmin, lonmax, wmin, wmax, itime_in);
            // cout << "tot_pwr: " << tot_pwr << "\n";

            // Or, try averaging first and integrating over a uniform patch:
            // (This matches the integration to 2 decimal places! At least for my single test. 11.17.16)            
            inp_pwr *= dlat * dlon * dw/8.;
            cout << "input power: " << inp_pwr << "\n";


            // Always do at least 2 steps in each axis (corner rays)
            num_freqs_fine = max(2, (int)floor( (wmax - wmin)/(2*PI*FREQ_STEP_SIZE )));
            num_lats_fine  = max(2, (int)floor( (dlat*1e-3)/(LAT_STEP_SIZE) ));
            num_lons_fine  = max(2, (int)floor( (dlon*1e-3)/(LON_STEP_SIZE) ));

            cout << "num steps: " << num_freqs_fine << ", " << num_lats_fine << ", " << num_lons_fine << "\n";

            crossing_log = fopen(crossingFileName.c_str(), "w");
            // FILE* area_log = fopen("/shared/users/asousa/WIPP/3dWIPP/area_log.txt","w");

            // cout << "T_STEP: " << TIME_STEP << "\n";
            double hit_counter = 0;
            double crossing_counter = 0;


            // --------------------- Interpolate + look for crossings ------------------
            //                            ( The main event)            
            // -------------------------------------------------------------------------
            cout << "checking for crossings...\n";
            time(&run_tstart);

            // Interpolate the first frames:
            for (int zz=0; zz<8; zz++) { interp_rayF(cur_rays[zz], &(prev_frames[zz]), 0); }

            // Get input area at top of ionosphere:
            initial_area = polygon_frame_area(prev_frames);


            // Step forward in time:
            for (double tt=TIME_STEP; tt < tmax; tt+=TIME_STEP) {

                // interpolate current frames:
                for (int zz=0; zz<8; zz++) { interp_rayF(cur_rays[zz], &(cur_frames[zz]), tt); }

                    // double fa = polygon_frame_area(cur_frames);
                    // printf("t= %g, Area= %2.2f\n",tt, fa);
                    // fprintf(area_log, "%g %g %g\n",tt, fa);

                // Check damping of each ray, and abort if they're all below a threshold:
                bool below_damping_thresh = false;
                for (int zz=0; zz<8; zz++) {
                    // cout << cur_frames[zz].damping << "\n"; 
                    if (cur_frames[zz].damping < DAMPING_THRESH) { below_damping_thresh = true;}
                }

                if (below_damping_thresh) { 
                    cout << "Below damping threshold! tt= " <<  tt << "\n";
                    break;
                }

                // Check each EA segment:
                // #pragma omp parallel for
                for (int rr = 0; rr < NUM_EA; rr++) {

                    // Ignore anything that looks way out of range
                    if (coarse_mask(cur_frames, prev_frames, EA_array[rr])) {
                        hit_counter ++;

                        // Calculate the geometric factor (spreading of guide rays)
                        frame_area = polygon_frame_area(cur_frames);
                        geometric_factor = initial_area/frame_area;  
                        cell_area = FREQ_STEP_SIZE*1.0/(num_lons_fine*num_lats_fine*num_freqs_fine);
                        

                        // Interpolate on fine-scale grid:
                        // #pragma omp parallel for
                        for (double ii=0; ii < 1; ii+=1./num_lons_fine) {         
                            for (double jj=0; jj < 1; jj+= 1./num_lats_fine) {     
                                for (double kk=0; kk < 1; kk+= 1./num_freqs_fine) { 
                                    
                                    // Clear previous values
                                    r_cur =  {};
                                    r_prev = {};

                                    // (to do: Save r_curs to avoid having to recalculate it)
                                    interp_ray_positions(cur_frames, ii, jj, kk, &r_cur);
                                    interp_ray_positions(prev_frames,ii, jj, kk, &r_prev);

                                    // Bam -- we finally have some little rays to check for crossings.
                                    if (crosses_EA(r_cur.pos, r_prev.pos, EA_array[rr])) {

                                        crossing_counter++;

                                        interp_ray_data(cur_frames, ii, jj, kk, &r_cur);
                                        interp_ray_data(prev_frames,ii, jj, kk, &r_prev);

                                        // Write crossing to log (for plotting)
                                        fprintf(crossing_log, "%g %g %g %g %g %g\n",
                                            r_cur.pos[0], r_cur.pos[1], r_cur.pos[2], 
                                            r_prev.pos[0], r_prev.pos[1], r_prev.pos[2]);

                                        // store time and frequency for the middle of this interpolation
                                        r_cur.dt = (r_cur.time - r_prev.time);
                                        r_cur.dlat = dlat;
                                        r_cur.dlon = dlon;
                                        // r_cur.ds   = (r_cur.pos - r_prev.pos).norm()*R_E;   // Jacob uses ds between the EA segments... hm

                                        t_grid = floor(r_cur.time/(TIME_STEP));
                                        f_grid = floor(kk*num_freqs_fine); 
                                                                           
                                        grid_ind = make_pair(t_grid, f_grid);

                                        cellT cell_cur = new_cell(r_cur);

                                        cell_cur.Lsh = EA_array[rr].Lsh;
                                        cell_cur.lat = EA_array[rr].lat;

                                        // Total power within this cell
                                        cell_cur.pwr = pow(inp_pwr * cell_area * geometric_factor * r_cur.damping, 2);

                                        if (crossing_db[rr].count(grid_ind)==0) {
                                            // If we haven't hit this same (time, freq, EA) combo yet, add it:
                                            crossing_db[rr].insert(make_pair(grid_ind, cell_cur));
                                        } else {
                                            // Else, sum the current frame with previous frames, so we can average:
                                            add_cell(&(crossing_db[rr].at(grid_ind)), &cell_cur);
                                        }
                                    }   // Crossings
                                }   // kk
                            }   // jj
                        }   // ii
                    }   // Coarse mask
                }   // EA array (single fieldline)

                // Step forward one frame:
                for (int zz=0; zz<8; zz++) { prev_frames[zz] = cur_frames[zz]; }

                } // tt

            cout << "hit counter: " << hit_counter << "\n";
            cout << "crossing counter: " << crossing_counter << "\n";


            }   // Distance from flash     

            time(&run_tend);
        } // Cur_rays

        cout << "crossing detection took " << (run_tend - run_tstart) << " sec\n";

        // Calculate scattering at crossings:
        for (int rr=NUM_EA-1; rr>=0; rr--) {
            calc_resonance(crossing_db[rr], EA_array[rr], da_N, da_S);
        }

        // // Step thru to the next frequency (shallow copy)
        // raylist_low = raylist_hi;

        // Write pN, pS files:
        cout << "Saving pN, pS files\n";

        ostringstream pN_name, pS_name;

        pN_name << out_dir << "/pN_" << out_lat << "_" << out_lon << "_" << round(wmin/(2*PI)) << ".dat";
        pS_name << out_dir << "/pS_" << out_lat << "_" << out_lon << "_" << round(wmin/(2*PI)) << ".dat";
        cout << pN_name.str() << "\n";
        cout << pS_name.str() << "\n";
        write_p_array(da_N, pN_name.str());
        write_p_array(da_S, pS_name.str());


    return 0; // Return statement.
} // Closing Main.