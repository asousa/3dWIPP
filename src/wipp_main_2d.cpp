#include <wipp.h>
using namespace std;
using namespace Eigen;

// #pragma omp 
int main(int argc, char *argv[]) 
{
    // map <int, rayF> raylist_hi;
    // map <int, rayF> raylist_low;
    map <int, rayF> raylist;

    double flash_pos[3];  // Flash position, geographic spherical coords (r, lat, lon)
    double flash_pos_sm[3];
    double x_in[3];
    double x_out[3];

    double flash_I0;
    double f1, f2;

    string ray_inp_dir;
    string out_dir;

    // debugging and dump filenames
    string dumpFileName;
    ostringstream eaFileName;
    string crossingFileName;

    // string low_file, high_file;

    int itime_in[2];
    
    double tmp_coords[3];
    double tmp_coords2[3];

    double in_lat, in_lon;
    double avg_distance_from_flash;

    double x[3];
    double lat0, lon0, rad0;
    // double maglat0, maglon0, magrad0;

    double dd;

    double wmin, wmax, latmin, latmax, lonmin, lonmax;
    double dw, dlon, dlat;  
    int tmax;

 
    rayT r_cur;    // Current interpolated ray (single timestep)
    rayT r_prev;   // Previous interpolated ray

    // Vector3d l0, l1;
    double* start_pos;

    // rayF* ray;
    rayF ray;

    // // // Output arrays (change in pitch angle, north and south hemispheres)
    // double da_N[NUM_E][NUM_TIMES] = {0};
    // double da_S[NUM_E][NUM_TIMES] = {0};

    // map <pair<int, int>, cellT> crossing_db[NUM_EA];
    // map <pair<int, int>, cellT> crossing_db;
    pair<int, int> grid_ind;
    int t_grid, f_grid;

    map <int, vector<double> > start_locs;
    // vector < vector<int> > adjacency_list;  // ( n x 4 ) list of adjacent ray indexes
    vector < vector<double> > adjacent_rays;
    // Array of pointers to current rays of interest
    // rayF* cur_rays[8];
    rayF cur_rays[8];
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
   
    vector<EA_segment> EA_array;

    // For calculating longitude variation -- pseudo-3d method
    vector<double> inp_pwrs;
    vector<double> offset_lons;
    

    // Location to determine output at (geomagnetic)
    double out_lat = 50;
    double out_lon = 0;
    double targ_lon = 0;

    int opt = 0;
    int opt_index = 0;

    double MLT_in = 0;
    double lon_spacing = 1;  // deg
    double num_lons = 1;    // number of longitudes on either side of the flash meridian
    // Default parameters:
    out_dir = "/shared/users/asousa/WIPP/3dWIPP/outputs/";
    dumpFileName= "fieldline_dump.dat";
    crossingFileName = "crossing_log.txt";
    flash_pos = {1, 0, 0};
    flash_I0 = -100e3;
    out_lat = 45;
    out_lon = 0;
    iyr = 2010; idoy = 1; isec = 0;
    // low_file  = "/shared/users/asousa/WIPP/3dWIPP/outputs/four_adjacent/rayout_1000_damped.ray";
    // high_file = "/shared/users/asousa/WIPP/3dWIPP/outputs/four_adjacent/rayout_1100_damped.ray";
    ray_inp_dir = "/shared/users/asousa/WIPP/3dWIPP/outputs/rays/";


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
        {"f1",          required_argument,    0, 'k'},
        {"f2" ,         required_argument,    0, 'l'},
        {"ray_dir",     required_argument,    0, 'm'},
        {"b_model",     required_argument,    0, 'n'},
        {"mlt",         required_argument,    0, 'o'},
        {"lon_spacing", required_argument,    0, 'p'},
        {"num_lons",    required_argument,    0, 'q'},
        {"num_freq_steps", required_argument, 0, 'r'},
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
                f1 = strtod(optarg, NULL);              break;
            case 'l':   // upper rayfile
                f2 = strtod(optarg, NULL);              break; 
            case 'm':
                ray_inp_dir = (string) optarg;          break;
            case 'n':
                model_number = atoi(optarg);            break;
            case 'o':
                MLT_in = strtod(optarg, NULL);          break;
            case 'p':
                lon_spacing = strtod(optarg, NULL);     break;
            case 'q':
                num_lons    = strtod(optarg, NULL);     break;
            case 'r':
                num_freqs_fine = atoi(optarg);          break;
            case '?':
                 printf("\nUnknown option: %s\n",opt);  break;
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
    cout << "B-field model: " << model_number << "\n";
    cout << "f1: " << f1 << endl;
    cout << "f2: " << f2 << endl;
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
    // tmp_coords = {rad0, lat0, lon0};
    // sphcar(flash_pos, tmp_coords);
    // Get flash position in SM coords:
    pol_to_cart_d_(&lat0, &lon0, &rad0, tmp_coords);
    mag_to_sm_d_(itime_in, tmp_coords, flash_pos_sm);

    // Vector of longitude offsets:
    if (num_lons > 0) {
        for (double ll=lon0; ll <= lon0 + num_lons*lon_spacing; ll+=lon_spacing) {
            offset_lons.push_back(ll);
        }
    }
    cout << "doing offset longitudes: ";
    print_vector(offset_lons);
// ---------------- Find all available rays --------------------------
    ostringstream tmp;
    tmp << ray_inp_dir << "/f_" << f1;
    vector <vector<double> > available_rays;

    get_available_rays(tmp.str(), &available_rays);

    cout << "found " << available_rays.size() << " available rays\n";

    adjacent_rays = find_adjacent_rays_2d(available_rays);
    cout << "found " << adjacent_rays.size() << " sets of adjacent rays\n";


    // Start with a fresh crossing db:
    // for (int rr=0; rr < EA_array.size(); ++rr) { crossing_db[rr].clear();} 
    map <pair<int, int>, cellT> crossing_db[EA_array.size()];


    // // -------------- Iterate through adjacent ray sets: -----------------

    for (vector < vector<double> >::iterator adj_row=adjacent_rays.begin(); adj_row!=adjacent_rays.end(); ++adj_row) {
        avg_distance_from_flash = 0;

        for (int i=0; i<2; ++i) {
            in_lat = (*adj_row)[2*i  ];
            in_lon = (*adj_row)[2*i+1];
            // cout << in_lat << "," << in_lon << endl;
            // in_lon = flash_pos[2]; // Just compare distance in latitude, since we'll rotate anyway
            avg_distance_from_flash += haversine_distance(in_lat, in_lon, flash_pos[1], flash_pos[2]);
        }
        avg_distance_from_flash /= 2000.0;  // Average dist of the 2 corner rays, in km
        
        if ( (avg_distance_from_flash <= MAX_GROUND_DISTANCE) ) {
            // print_vector(*adj_row);

            // -------------- Load current rays ----------------------------------:
            for (int jj=0; jj<8; ++jj) {
                ostringstream cur_file;
                double jj2 = fmod(1.0*jj, 2.);
                // f0, lat0, lon0
                // f0, lat1, lon0
                // f0, lat0, lon1
                // f0, lat1, lon1
                // f1, lat0, lon0
                // f1, lat1, lon0
                // f1, lat0, lon1
                // f1, lat1, lon1
                
                double f   = (jj < 4 ? f1 : f2);
                double lat = (jj2 < 2 ? (*adj_row)[2*jj2]   : (*adj_row)[2*(jj2 - 2)]);
                double lon = (jj2 < 2 ? (*adj_row)[2*jj2+1] : (*adj_row)[2*(jj2 - 2) + 1]);

                // cout << "jj: " << jj << ", jj2: " << jj2 << ", " << f << ", " << lat << ", " << lon << endl;
                cur_file << ray_inp_dir   << "/f_" << f   << "/lon_"    << lon
                         << "/ray_" << f  << "_"   << lat << "_" << lon << ".ray";

                raylist = read_rayfile(cur_file.str());
                cur_rays[jj] = raylist.at(1); // First ray should have ray number "1"
         
                calc_stix_parameters(&(cur_rays[jj]));

                // Also load the damping file: 
                cur_file.str(""); cur_file.clear();
                cur_file << ray_inp_dir   << "/f_" << f   << "/lon_"    << lon
                         << "/damp_" << f  << "_"   << lat << "_" << lon << ".ray";
                raylist = read_dampfile(cur_file.str());
                
                cur_rays[jj].damping = raylist.at(1).damping;

                // Rotate the center ray to the guide longitudes on either side:
                targ_lon = (fmod(jj, 4.) < 2 ? flash_pos[2] -1.*lon_spacing/2. : flash_pos[2] + lon_spacing/2. );
                for (int i=0; i < cur_rays[jj].pos.size(); ++i) {
                    double tmp_mag[3], tmp_sm[3];

                    sm_to_mag_d_(itime_in, cur_rays[jj].pos[i].data(),tmp_mag);
                    cardeg(tmp_mag);
                    tmp_mag[2] = targ_lon;
                    degcar(tmp_mag);
                    mag_to_sm_d_(itime_in, tmp_mag, cur_rays[jj].pos[i].data());
                }
                cur_rays[jj].in_lat = lat; cur_rays[jj].in_lon = targ_lon;

                cout << "cur_rays[" << jj <<"]: " << cur_rays[jj].w/2./PI << 
                      ", " << cur_rays[jj].in_lat << ", " << cur_rays[jj].in_lon << endl;
            }


            
            if (DEBUG) {check_memory_usage();}
        
            // Find minimum and maximum frequencies, start lats, and start lons:
            wmin   = cur_rays[0].w;         wmax = cur_rays[0].w;
            latmin = cur_rays[0].in_lat;  latmax = cur_rays[0].in_lat;
            lonmin = cur_rays[0].in_lon;  lonmax = cur_rays[0].in_lon;
            tmax   = cur_rays[0].time.back();
            double in_lat, in_lon;

            for (int i=1; i < 8; i++) {
                in_lon = cur_rays[i].in_lon;
                in_lat = cur_rays[i].in_lat;
                
                if (in_lon >= 360)                  { in_lon -= 360; }
                if (cur_rays[i].w < wmin)          { wmin = cur_rays[i].w; } 
                if (cur_rays[i].w > wmax )         { wmax = cur_rays[i].w; }
                if (cur_rays[i].in_lat < latmin )  { latmin = cur_rays[i].in_lat; }
                if (cur_rays[i].in_lat > latmax )  { latmax = cur_rays[i].in_lat; }
                if (in_lon < lonmin )               { lonmin = in_lon; }
                if (in_lon > lonmax )               { lonmax = in_lon; }
                if (cur_rays[i].time.back() < tmax){ tmax = cur_rays[i].time.back(); }
            }

            // starting separation in lat, lon directions (meters)
            dlat = D2R*(R_E + H_IONO)*(latmax - latmin);
            dlon = D2R*(R_E + H_IONO)*(lonmax - lonmin)*cos(D2R*(latmax + latmin)/2.);
            dw   = wmax - wmin;

            cout << "\n----- current rays: -----\n";
            cout << "lon: " << lonmax << ", " << lonmin << " deg\n";
            cout << "lat: " << latmax << ", " << latmin << " deg\n";
            cout << "f: " << wmax/(2*PI) << ", " << wmin/(2*PI) << " Hz\n"; 
            cout << "Avg distance from flash: " << avg_distance_from_flash << " km\n";


            // Scale the input power by dlat, dlon, dw:
            // (ray spacing may not be consistent)
            double inp_pwr = 0;
            inp_pwrs.clear();
            for (vector<double>::iterator offlon=offset_lons.begin(); offlon!=offset_lons.end(); ++offlon) {
                double tmp_lonmin = lonmin + *offlon;
                double tmp_lonmax = lonmax + *offlon;
                cout << tmp_lonmin << " " << tmp_lonmax << " " << *offlon << endl;
                inp_pwr = total_input_power(flash_pos_sm, flash_I0, 
                    latmin, latmax, tmp_lonmin, tmp_lonmax, wmin, wmax, itime_in);

                inp_pwrs.push_back(inp_pwr);

            }

            cout << "input power vec: ";
            print_vector(inp_pwrs);

            inp_pwr = total_input_power(flash_pos_sm, flash_I0, 
                                        latmin, latmax, lonmin, lonmax, wmin, wmax, itime_in);
            cout << "input energy: " << inp_pwr << " Joules\n";

            if (CROSSING_METHOD==0) {
            // Always do at least 1 step in each axis:
                num_freqs_fine = max(1, (int)floor( (wmax - wmin)/(2.*PI*FREQ_STEP_SIZE )));
                num_lats_fine  = max(1, (int)floor( (dlat*1e-3)/(LAT_STEP_SIZE) )); 
                num_lons_fine  = max(1, (int)floor( (dlon*1e-3)/(LON_STEP_SIZE) ));
                cout << "num steps: " << num_freqs_fine << ", " << num_lats_fine << ", " << num_lons_fine << "\n";
            }
            
            crossing_log = fopen(crossingFileName.c_str(), "w");

            // FILE* area_log = fopen("/shared/users/asousa/WIPP/3dWIPP/area_log.txt","w");

            double hit_counter = 0;
            double crossing_counter = 0;

          
            // --------------------- Interpolate + look for crossings ------------------
            //                            ( The main event)            
            // -------------------------------------------------------------------------
            cout << "checking for crossings...\n";
            time(&run_tstart);
  
            // Interpolate the first frames:
            for (int zz=0; zz<8; zz++) { interp_rayF(&cur_rays[zz], &(prev_frames[zz]), 0); }

            // Get input area at top of ionosphere:
            // initial_area = polygon_frame_area(prev_frames);


            // Step forward in time:
            for (double tt=TIME_STEP; tt < tmax; tt+=TIME_STEP) {
                // cout << "t: " << tt << "\n";
                // interpolate current frames:
                for (int zz=0; zz<8; zz++) { interp_rayF(&cur_rays[zz], &(cur_frames[zz]), tt); }

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
                // for (int rr = 0; rr < NUM_EA; rr++) {
                for (int rr = 0; rr < EA_array.size(); rr++) {

                    // Ignore anything that looks way out of range
                    if (coarse_mask(cur_frames, prev_frames, EA_array[rr])) {
                        hit_counter ++;
                        // cout << "hit" << endl;
                        // Calculate the geometric factor (spreading of guide rays)
                        frame_area = polygon_frame_area(cur_frames);
                        // geometric_factor = initial_area/frame_area;
                        // cell_area = FREQ_STEP_SIZE*1.0/(num_lons_fine*num_lats_fine*num_freqs_fine);

                        if (CROSSING_METHOD==1) {  
                            // Calculate crossings using the optmization method
                            // Vector2d soln;
                            // soln << 0.5, 0.5;
                            double ii_up, jj_up;
                            double ii_dn, jj_dn;
                            double ii, jj;
                            double S1, S2;
                            // find_crossing(cur_frames, prev_frames, EA_array[rr].ea_pos, 0, &ii_dn, &jj_dn);
                            // find_crossing(cur_frames, prev_frames, EA_array[rr].ea_pos, 1, &ii_up, &jj_up);
                            // cout << "lower f: " << ii_dn << ", " << jj_dn << "\n";
                            // cout << "upper f: " << ii_up << ", " << jj_up << "\n";

                            // if ((S1 >= 0 ) && (S1 <= 1) && (S2 >= 0) && (S2 <= 1)) {

                                // if ( (ii_dn >= 0 ) && (ii_dn <= 1) && (jj_dn >=0) && (jj_dn <=1) &&
                                     // (ii_up >= 0 ) && (ii_up <= 1) && (jj_up >=0) && (jj_up <=1) ) {
                                // cout << "lower f: " << ii_dn << ", " << jj_dn << "\n";
                                // cout << "upper f: " << ii_up << ", " << jj_up << "\n";

                                    for (double kk=0; kk < 1; kk += 1./num_freqs_fine) {

                                        find_crossing(cur_frames, prev_frames, EA_array[rr].ea_pos, kk, &ii, &jj);
                                        // ii = (1.-kk)*ii_dn + kk*ii_up;
                                        // jj = (1.-kk)*jj_dn + kk*jj_up;
                                        // if ( (ii >= 0 ) && (ii <= 1) && (jj >=0) && (jj <=1) && (S1 >= 0) && (S2 <= 1)) {
                                        if ( (ii >= 0 ) && (ii <= 1) && (jj >=0) && (jj <=1)) {
                                            r_cur =  {};
                                            r_prev = {};
                                            
                                            interp_ray_positions(cur_frames, ii, jj, kk, &r_cur);
                                            interp_ray_positions(prev_frames,ii, jj, kk, &r_prev);
                                            interp_ray_data(cur_frames, ii, jj, kk, &r_cur);
                                            interp_ray_data(prev_frames,ii, jj, kk, &r_prev);
                                            

                                            // // Write crossing to log (for plotting)
                                            // crossing_log = fopen("crossing_log_newway.txt", "a");
                                            // fprintf(crossing_log, "%g %g %g %g %g %g\n",
                                            //     r_cur.pos[0], r_cur.pos[1], r_cur.pos[2], 
                                            //     r_prev.pos[0], r_prev.pos[1], r_prev.pos[2]);
                                            // fclose(crossing_log);

                                            t_grid = floor(tt/TIME_STEP); //floor(r_cur.time/(TIME_STEP));
                                            f_grid = floor(r_cur.w); //floor(kk*num_freqs_fine); 
      
                                            grid_ind = make_pair(t_grid, f_grid);
                                            cellT cell_cur = new_cell(r_cur);

                                            // Get length between planes for volume calc:
                                            double frame_length = (r_cur.pos - r_prev.pos).norm()*R_E;
                                            // Volume of cell, in meters
                                            double cell_vol = frame_area*frame_length;

                                            // cell_cur.pwr = (inp_pwr / frame_area)*(1.0*FREQ_STEP_SIZE/num_freqs_fine)*(r_cur.damping);
                                            cell_cur.pwr = (inp_pwr / cell_vol)*(1.0/num_freqs_fine)*(r_cur.damping);

                                            // current power for adjacent longitude bins

                                            for (vector<double>::iterator inpwr=inp_pwrs.begin(); 
                                                                inpwr!=inp_pwrs.end(); ++inpwr) {
                                                cell_cur.pwr_vec.push_back((*inpwr / frame_area)*
                                                    (1.0*FREQ_STEP_SIZE/num_freqs_fine)*(r_cur.damping));
                                            }

                                            // cout << "t: " << t_grid << " f: " << f_grid/2./PI;
                                            // cout << " energy density vector: ";
                                            print_vector(cell_cur.pwr_vec);
                                            // cout << " cell energy density: " << cell_cur.pwr << " J/m^3\n";
                                            crossing_db[rr].insert(make_pair(grid_ind, cell_cur));
                                    }
                                }
                            // }
                        } else {
                            // Calculate crossings by interpolating and checking

                            // Interpolate on fine-scale grid:
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
                                            crossing_log = fopen("crossing_log_oldway.txt", "a");
                                            fprintf(crossing_log, "%g %g %g %g %g %g\n",
                                                r_cur.pos[0], r_cur.pos[1], r_cur.pos[2], 
                                                r_prev.pos[0], r_prev.pos[1], r_prev.pos[2]);
                                            fclose(crossing_log);

                                            // store time and frequency for the middle of this interpolation
                                            // r_cur.dt = (r_cur.time - r_prev.time);
                                            r_cur.dlat = dlat;
                                            r_cur.dlon = dlon;
                                            // r_cur.ds   = (r_cur.pos - r_prev.pos).norm()*R_E;   // Jacob uses ds between the EA segments... hm

                                            // t_grid = floor(r_cur.time/(TIME_STEP));
                                            // f_grid = floor(kk*num_freqs_fine); 

                                            t_grid = floor(tt/TIME_STEP); //floor(r_cur.time/(TIME_STEP));
                                            f_grid = floor(r_cur.w); //floor(kk*num_freqs_fine); 

                                                                                   
                                            grid_ind = make_pair(t_grid, f_grid);
                                            // cout << "grid inds: " << t_grid << ", " << f_grid << "\n";
                                            cellT cell_cur = new_cell(r_cur);

                                            // cell_cur.Lsh = EA_array[rr].Lsh;
                                            // cell_cur.lat = EA_array[rr].lat;

                                            // Total power within this cell
                                            // cell_cur.pwr = pow(inp_pwr * geometric_factor * cell_area
                                            //      * r_cur.damping, 2);

                                            // (total power / frame area)
                                            //      * (damping losses @ this cell)
                                            //      * (cell size in frequency axis)
                                            // cout << "inp pwr: " << inp_pwr;
                                            // cout << " frame area: " << frame_area;
                                            // cout << " other factor: " << (1.0*FREQ_STEP_SIZE/num_freqs_fine);
                                            // cout << " damping: " << r_cur.damping;

                                            cell_cur.pwr = (inp_pwr / frame_area)*(1.0*FREQ_STEP_SIZE/num_freqs_fine)*(r_cur.damping);
                                            // cout << "t: " << t_grid << " f: " << f_grid;
                                            // cout << " cell pwr: " << cell_cur.pwr << "\n";

                                            // // current power for adjacent longitude bins
                                            // for (vector<double>::iterator inpwr=inp_pwrs.begin(); 
                                            //                     inpwr!=inp_pwrs.end(); ++inpwr) {
                                            //     cell_cur.pwr_vec.push_back((*inpwr / frame_area)*
                                            //         (1.0*FREQ_STEP_SIZE/num_freqs_fine)*(r_cur.damping));
                                            // }

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
                            }   // ii */
                        } // Detection mode
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

    // Loop over offset longitudes
    for (int offset_index=0; offset_index < offset_lons.size(); ++offset_index) {
        cout << "doing offset " << offset_index << ", offset_lon: " << offset_lons[offset_index] << endl;
        // Clear output arrays
        // // Output arrays (change in pitch angle, north and south hemispheres)
        double da_N[NUM_E][NUM_TIMES] = {0};
        double da_S[NUM_E][NUM_TIMES] = {0};

        // Calc resonance at each EA segment:
        // for (int rr=NUM_EA-1; rr>=0; rr--) {
        for (int rr=EA_array.size()-1; rr>=0; rr--) {
            calc_resonance(crossing_db[rr], EA_array[rr], da_N, da_S, offset_index);
        }

        // Write pN, pS files:
        cout << "Saving pN, pS files\n";

        ostringstream pN_name, pS_name;

        pN_name << out_dir << "/pN_" << out_lat << "_" << offset_lons[offset_index] << "_" << round(wmin/(2*PI)) << ".dat";
        pS_name << out_dir << "/pS_" << out_lat << "_" << offset_lons[offset_index] << "_" << round(wmin/(2*PI)) << ".dat";
        cout << pN_name.str() << "\n";
        cout << pS_name.str() << "\n";
        write_p_array(da_N, pN_name.str());
        write_p_array(da_S, pS_name.str());
    }



    return 0; // Return statement.
} // Closing Main.