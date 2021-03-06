#include <wipp.h>

using namespace std;
using namespace Eigen;


// double input_power_scaling(double* flash_loc, double* ray_loc, double mag_lat, double w, double i0) {
//     // Returns ray power at the top of the ionosphere
//     // per unit in area and frequency.

//     double theta;        // angle between two vectors
//     double gc_distance; // great circle distance between two points
//     double dist_tot;
//     double xi;

//     double S, S_vert;
//     double attn_factor;
//     double w_sq, f;

//     // Vector3d v1;
//     // Vector3d v2;
//     double v1[3];
//     double v2[3];

//     f = w/(2.0*PI);
    
//     cardeg(flash_loc, v1);
//     cardeg(ray_loc,   v2);
//     gc_distance = haversine_distance(v1[1], v1[2], v2[1], v2[2]);
//     // v1 = Map<VectorXd>(flash_loc,3,1);
//     // v2 = Map<VectorXd>(ray_loc,  3,1);

//     // theta = acos(v1.dot(v2)/(v1.norm()*v2.norm()));
//     // // Arc length (~great-circle distance) between vI0ectors
//     // gc_distance = (R_E)*theta;

//     // cout << "gc_distance: " << gc_distance << "\n";
//     // total distance up to ionosphere:
//     dist_tot = hypot(gc_distance, H_IONO);
//     xi = atan2(gc_distance, H_IONO);  // Incident angle

//     w_sq =  pow( w , 2 );
//     S = ( (1/Z0) * pow( (H_E*i0*2E-7*(sin(xi)/dist_tot)*w*(P_A-P_B)) , 2 ) 
//                    /  (  (w_sq+pow(P_A,2))*(w_sq+pow(P_B,2))  )      ) ;
//     S_vert = S * cos(xi) ;  // factor for vert prop.

//     // Ionosphere absorption model
//     attn_factor = pow(10,-(ionoAbsorp(mag_lat,f)/10)  );
//     S_vert = S_vert * attn_factor;

//     // printf("i0: %2.3f, dist_tot: %2.3f, xi: %2.3f, S_vert; %e\n",i0, dist_tot, xi, S_vert);
    
//     return S_vert;
// }




// double total_input_power(double flash_pos_sm[3], double i0, 
//                         double latmin, double latmax, double lonmin, double lonmax, double wmin, double wmax, int itime_in[2]) {
//     // Determine the total input power tracked by the set of guide rays:

//     double tot_pwr = 0;
//     double pwr = 0;
//     // Integration step sizes
//     double dlat = 0.05; 
//     double dlon = 0.05;
//     double dw   = 5*2*PI;
//     double tmp_coords[3] = {0,0,0};
//     double x_sm[3];

//     for (double w = wmin + dw/2; w < wmax; w+=dw) {
//         for (double lat = latmin + dlat/2; lat < latmax; lat+=dlat) {
//             for (double lon=lonmin; lon < lonmax; lon+=dlon) {
//                 // cout << "(" << lat << ", " << lon << ")\n";
//                 tmp_coords = {1 + H_IONO/R_E, lat, lon};
//                 degcar(tmp_coords);
//                 mag_to_sm_d_(itime_in, tmp_coords, x_sm);

//                 pwr = input_power_scaling(flash_pos_sm, x_sm, lat, w, i0);

//                 double dist_lat = (R_E + H_IONO)*dlat*D2R;
//                 double dist_lon = (R_E + H_IONO)*dlon*sin(D2R*lat)*D2R;
//                                 // Latitude distance      longitude distance       freq dist
//                 // cout << "dist_lat: " << dist_lat << ", dist_lon: " << dist_lon << "\n";
//                 tot_pwr += pwr * dist_lat * dist_lon * dw;
//             }
//         }
//     }

//     return tot_pwr;
// }









void interp_ray_positions(rayT framelist[8],  double n_x, double n_y, double n_z, rayT* rayout) {
    // Fine-scale interpolation (positions only)

    double W[8];
    W[0] = (1. - n_x)*(1. - n_y)*(1. - n_z);
    W[1] = n_x*(1. - n_y)*(1. - n_z);
    W[2] = n_y*(1. - n_x)*(1. - n_z);
    W[3] = n_x*n_y*1.*(1. - n_z);
    W[4] = (1. - n_x)*(1. - n_y)*n_z;
    W[5] = n_x*(1. - n_y)*n_z;
    W[6] = n_y*(1. - n_x)*n_z;
    W[7] = n_x*n_y*n_z*1.;

    rayout->pos = {0, 0, 0};

    for (int jj=0; jj<8; jj++){  // Corner rays
        for (int ii=0; ii<3; ii++){  // X, Y, Z
            (rayout->pos)[ii]   += W[jj]*((framelist[jj].pos).data()[ii]);
        }
    }

    rayout->time = framelist[0].time;
}



// void interp_ray_data(rayF** raylist, double n_x, double n_y, double n_z, int t_ind, rayT* rayout) {
void interp_ray_data(rayT framelist[8], double n_x, double n_y, double n_z, rayT* rayout) {

    // Interpolate the rest of the stuff we'll need in the ray

    double W[8];
    W[0] = (1. - n_x)*(1. - n_y)*(1. - n_z);
    W[1] = n_x*(1. - n_y)*(1. - n_z);
    W[2] = n_y*(1. - n_x)*(1. - n_z);
    W[3] = n_x*n_y*1.*(1. - n_z);
    W[4] = (1. - n_x)*(1. - n_y)*n_z;
    W[5] = n_x*(1. - n_y)*n_z;
    W[6] = n_y*(1. - n_x)*n_z;
    W[7] = n_x*n_y*n_z*1.;

    for (int jj=0; jj<8; jj++){  // Corner rays

        // Vector-valued
        for (int ii=0; ii<3; ii++){  // X, Y, Z
            // (rayout->pos)[ii]   += W[jj]*((raylist[jj]->pos[t_ind]).data()[ii]);
            // (rayout->vgrel)[ii] += W[jj]*((raylist[jj]->pos[t_ind]).data()[ii]);
            (rayout->n)[ii]     += W[jj]*((framelist[jj].n).data()[ii]);
            (rayout->B0)[ii]    += W[jj]*((framelist[jj].B0).data()[ii]);

        }

        // scalar-valued here

        rayout->w += W[jj]*(framelist[jj].w);

        rayout->inp_pwr += W[jj]*(framelist[jj].inp_pwr);
        rayout->damping += W[jj]*(framelist[jj].damping);

        rayout->stixR += W[jj]*(framelist[jj].stixR);
        rayout->stixL += W[jj]*(framelist[jj].stixL);
        rayout->stixP += W[jj]*(framelist[jj].stixP);

        rayout->in_lat+= W[jj]*(framelist[jj].in_lat);
        rayout->in_lon+= W[jj]*(framelist[jj].in_lon);
    }
}


void calc_stix_parameters(rayF* ray) {
    // Calculate Stix parameters for an entire ray (all timesteps)
    // Results are stored in the ray structure.

    double w, whs, wps2;
    double R, L, P, S, D, A, B;
    double kpar, kperp, kmag;
    double theta;
    double B0mag;
    double sin_th, cos_th, sin_th_sq, cos_th_sq;


    Vector3d B0;
    Vector3d Bhat;
    Vector3d n_vec;
    Vector3d k;
    // cout << "Stix params... \n";
    for (int ii=0; ii < ray->time.size(); ii++) {
        // ---------- Evaluate Stix parameters: ------
        // cout << "t: " << ii << "\n"; 
        wps2 = 0;
        R = 1.;
        L = 1.;
        P = 1.;

        w = ray->w;
        B0    = Map<VectorXd>(ray->B0[ii].data(), 3,1);
        n_vec = Map<VectorXd>(ray->n[ii].data(),3,1);

        B0mag = B0.norm();


        // Sum over constituents
        for (int jj=0; jj < ray->Ns[ii].size(); jj++) {
            
            // Ns*(Qs^2)/(ms*Eps_0)
            wps2 = ray->Ns[ii][jj]*pow(ray->qs[jj],2) \
                   /(ray->ms[jj]*EPS0);
            // qB/m
            whs  = ray->qs[jj]*B0mag/ray->ms[jj];

            // Complex modification to whs -- for now, ignore. (8.19.2016)
            // wcs  = whs * ray.w/(ray.w + 1i*ray.nus[ii][jj]);

            R-= wps2/(w*(w + whs));
            L-= wps2/(w*(w - whs));
            P-= wps2/(w*w);
        }
        S = (R + L)/2.;
        D = (R - L)/2.;

        ray->stixR.push_back(R);
        ray->stixL.push_back(L);
        ray->stixP.push_back(P);
        // ray->stixS.push_back(S);
        // ray->stixD.push_back(D); 
        // --------------------------------------------

    } // ii (timestep)
}


vector<EA_segment> init_EA_array(double lat, double lon, int itime_in[2], int model_number) {

    // vector<EA_segment> EA_array(NUM_EA);
    vector<EA_segment> EA_array;

    double x_in[3], x_in_magcar[3], x_sm[3];         
    double x_cur[3], x_prev[3], x_out[3];

    double w1, w2;  // Interpolation weights
    double prev_lat;

    double Bo[3];
    double Bomag;

    double x_fl[TRACER_MAX][3]; // elements along the field line
    double b_fl[TRACER_MAX];    // b-field magnitude along line
    double b_dipole[3], b_sm[3];    
    int Nsteps;

    double lam, slam, clam, dL_lam, ptR, ptX, ptY, x1, x2, y1, y2;
    double slam2, clam2, rootTerm, x_unit_vect, y_unit_vect, ptL;
    double slat_term;
    double EA_a, EA_b, EA_c;
    double x_iono[3], x_iono_tmp[3];
    double B_iono[3], B_eq[3];
    double alpha_eq;

    double dx;
    double tmp_ftc;
    int ind;



    x_in = {1, lat, lon}; // Start at the Earth; we'll select bins above H_MAGNETO down below.

    // Get start coordinate in SM cartesian:
    x_in_magcar = {x_in[0], x_in[1], x_in[2]};
    degcar(x_in_magcar);
    mag_to_sm_d_(itime_in, x_in_magcar, x_sm);

    // cout << "\torig: ";
    // print_array(x_in, 3);
    // cout << "\tSM: ";
    // print_array(x_sm, 3);

    
    // This works!
    //Nsteps = trace_fieldline(x_in, x_fl, TRACER_STEP);

    double tsyg_params[10] = {0};
    double VG[3];
   
    // Setup for IGRF:    
    load_TS05_params(itime_in, tsyg_params, VG);
    init_igrf(itime_in);

    // set DST:
    tsyg_params[1] = -20;

    // Trace field line:
    Nsteps = trace_fieldline(itime_in, x_sm, x_fl, b_fl, TRACER_STEP, model_number, tsyg_params);


    double lats[Nsteps];
    double lons[Nsteps];
    double dist_n[Nsteps];
    double ftc_n[Nsteps];   // Flight time constant to northern hemi

    int EA_index = 0;
    
    // Get magnetic latitudes of entries
    for (int i=0; i < Nsteps; i++) {
        sm_to_mag_d_(itime_in, x_fl[i], x_cur);
        cardeg(x_cur);
        lats[i] = x_cur[1];
        lons[i] = x_cur[2];
        // Store the cumulative distance to northern hemisphere:
        if (i==0) {
            dist_n[i] = 0;
        } else {
            dist_n[i] = sqrt(pow(x_fl[i][0] - x_fl[i-1][0], 2)
                            +pow(x_fl[i][1] - x_fl[i-1][1], 2)
                            +pow(x_fl[i][2] - x_fl[i-1][2], 2))
                        + dist_n[i - 1];
        }
        // cout << "d[" << i << "] " << dist_n[i] << "\n";

    }

    // cout << "total distance: " << dist_n[Nsteps] << "\n";
    // Get effective L-shell
    // (Do we take this to be the radius at geomag equator, or the maximum?)
    int Lsh_index = nearest(lats, Nsteps, 0, true);

    double Lsh = norm(x_fl[Lsh_index], 3);
    cout << "\tL shell: " << Lsh << "\n";

    // Get b-field at equator:
    bmodel(itime_in, x_fl[Lsh_index], tsyg_params, model_number, B_eq);

    // Get b-field at H_IONO:
    sm_to_mag_d_(itime_in, x_fl[0], x_iono_tmp);
    cardeg(x_iono_tmp);
    x_iono_tmp[0] += H_IONO/R_E;
    degcar(x_iono_tmp);
    mag_to_sm_d_(itime_in, x_iono_tmp, x_iono);
    bmodel(itime_in, x_iono, tsyg_params, model_number, B_iono);

    // get loss cone at equator:
    alpha_eq = asin(sqrt(norm(B_eq,3)/norm(B_iono, 3)));
    // cout << "Alpha_eq: " << alpha_eq << "\n";

    // Calculate flight-time constants:
    int iono_ind = 0;
    double alt = 0;
    // Find index of ionosphere
    while (alt < (H_IONO/R_E + 1)) {
        iono_ind +=1;
        Vector3d xtmp = Map<VectorXd>(x_fl[iono_ind],3,1);
        alt = xtmp.norm();
    }

    // cout << "Iono ind: " << iono_ind << "\n";

    // Integrate Walt 4.25
    for (int i=0; i < Nsteps; i++) {
        if (i <= iono_ind) {
            ftc_n[i] = 0;
        } else {
            dx = R_E*(dist_n[i] - dist_n[i-1]);
            tmp_ftc = dx/(sqrt(1 - b_fl[i]/b_fl[iono_ind]) );
            // Skip over infs and nans (only occurs at b_fl = b_iono)
            if (isnan(tmp_ftc)) {
                ftc_n[i] = ftc_n[i-1];
            } else {
                ftc_n[i] = ftc_n[i-1] + tmp_ftc;    
            }
        }
        // cout << "i: " << i << " lat: " << lats[i] << " dist: " << dx <<" ftc: " << ftc_n[i] << "\n";
    }

    // // Evaluate Walt's flight-time constant, from equator to mirror pt.
    // double walt_tau = (0.117/2.)*Lsh*C*(1 - 0.4635*pow(sin(alpha_eq),0.75));

    // cout << "Walt tau: " << walt_tau << "\n";
    // cout << "Integrated tau: " << ftc_n[Nsteps-1] << "\n";



    // double targ_lat = EALimN;
    // Starting latitude
    double start_lat = acos(sqrt((1. + H_MAGNETO/R_E)/Lsh))*R2D;
    double stop_lat  =-1.0*start_lat;  
    start_lat = floor(start_lat/EAIncr)*EAIncr;
    double targ_lat = start_lat;
//  generate entries for each EA segment

    // for (int i=0; i < NUM_EA; i++) {
    while (targ_lat > (stop_lat)) {
        ind = nearest(lats, Nsteps, targ_lat, true);

        EA_segment seg;

        // cout << "lat: " << lats[ind] << "\n";
        // cout << "Lsh: " << Lsh << "\n";
        // L shell:
        seg.Lsh = Lsh;

        // Latitude:
        seg.lat = lats[ind];
        seg.lon = lons[ind];

        // Distance to northern and southern ionosphere intersections:
        seg.dist_to_n = dist_n[ind]                    - H_IONO/R_E;
        seg.dist_to_s = dist_n[Nsteps-1] - dist_n[ind] - H_IONO/R_E;

        // Flight-time constants:
        seg.ftc_n     = ftc_n[ind];
        seg.ftc_s     = ftc_n[Nsteps-1] - ftc_n[ind];
        // cout << "ftc n: " << EA_array[i].ftc_n << " ftc s: " << EA_array[i].ftc_s << "\n";

        // Intersection point:
        seg.ea_pos = Map<VectorXd>(x_fl[ind],3,1);
        // Get B:
        // Bomag = b_fl[ind];
        // cout << "Bo: " << Bomag << "\n";
        bmodel(itime_in, x_fl[ind], tsyg_params, model_number, Bo);

        // Get unit vector pointing along field line:
        // (i.e., normal vector to this EA segment)
        Bomag = norm(Bo, 3);
        // Bomag = b_fl[ind];
        seg.ea_norm = Map<VectorXd>(Bo, 3, 1)/Bomag;

        // Calculate the width in the same way Jacob did:
        // (Width in L_MARGIN, assuming a dipole field)
        lam = lats[ind];

        clam = cos(lam*D2R);
        slam = sin(lam*D2R);
        clam2 = pow(clam,2);
        slam2 = pow(slam,2);
        rootTerm = sqrt(1+3*slam2);
          
        dL_lam = clam2*clam / rootTerm * L_MARGIN ;
        
        x_unit_vect = (3*clam2 - 2) / rootTerm ;
        y_unit_vect = (3*slam*clam) / rootTerm ;
        
        ptR = ptL*clam2;
        ptX = ptR*clam;
        ptY = ptR*slam;
          
        x1 = ptX - x_unit_vect*dL_lam ;
        x2 = ptX + x_unit_vect*dL_lam ;
        
        y1 = ptY - y_unit_vect*dL_lam ;
        y2 = ptY + y_unit_vect*dL_lam ;
        
        EA_a = y1 - y2;
        EA_b = x2 - x1;
        EA_c = x1*y2 - y1*x2;
        
        seg.radius = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))/2.0; 
        // printf("EA_length: %g\n",2*R_E*EA_array[i].radius);


        // Calculate electron gyrofrequency:
        seg.wh = (Q_EL/M_EL)*Bomag;

        // Spatial derivative of gyrofrequency:
        // (Currently evaluating using Jacob's formula -- dipole model assumption!!)
        slat_term = sqrt(1+3*slam*slam);

        seg.dwh_ds = 3*(seg.wh)/(Lsh*R_E)*slam/slat_term*
                (1/(slat_term*slat_term) + 2/(clam*clam));



        // Calculate local loss cone:
        // double alpha_lc_calc = asin(sqrt( slat_term/pow(clam,6) )*sin(alpha_eq));
        // EA_array[i].alpha_lc = asin(sqrt( slat_term/pow(clam,6) )*sin(alpha_eq));
        seg.alpha_lc = asin(sqrt(Bomag/norm(B_iono,3)));
        
        // printf("lat: %g alpha_lc (mine): %g alpha_lc (theirs): %g\n",
        //     lats[ind], R2D*EA_array[i].alpha_lc, R2D*alpha_lc_calc);
        


        // Save equatorial loss cone:
        seg.alpha_eq = alpha_eq;

        // dv_|| / ds:  (again, Jacob's formula -- dipole model assumption)
        seg.dv_para_ds = -0.5*pow(sin(seg.alpha_lc),2)/
                                      cos(seg.alpha_lc)/
                                      seg.wh*seg.dwh_ds;

        // Ratio of B-field at equator vs local: (used in mapping pitch angles to equator)
        seg.Bo_ratio = sqrt(norm(B_eq,3)/norm(Bo,3));

        // Area of EA disc:
        seg.area = pow(R_E*seg.radius,2)*PI;


        // cout << "ratio: " << EA_array[i].Bo_ratio << "\n";
        // // analytical (for comparison)
        // double epsm = (1/Lsh)*(R_E+H_IONO)/R_E;
        // double alpha_eq = asin(sqrt( pow(epsm,3)/sqrt(1+3*(1-epsm)) ));
        // double slat_term = sqrt(1+3*slam*slam);
        // double alpha_calc = asin(sqrt( slat_term/pow(clam,6) )*sin(alpha_eq));

        // cout << "alpha_lc: " << EA_array[i].alpha_lc << " calculated: " << alpha_calc << "\n";

        EA_array.push_back(seg);
        // Bump index
        targ_lat-= EAIncr;
    }



    // Loop through again to calculate ds:
    // (Distance along field line between EA segments)
    for (int i=0; i < EA_array.size(); i++) {
        if (i==0) {
            EA_array[i].ds = EA_array[i+1].dist_to_n - EA_array[i].dist_to_n;
        } else if (i== EA_array.size() - 1) {
            EA_array[i].ds = EA_array[i].dist_to_n - EA_array[i-1].dist_to_n;
        } else {
            EA_array[i].ds = (EA_array[i+1].dist_to_n - EA_array[i].dist_to_n)*0.5
                            +(EA_array[i].dist_to_n - EA_array[i-1].dist_to_n)*0.5;
        }

        EA_array[i].ds *= R_E;  // Meters
    }



    return EA_array;
}







void dump_fieldlines(int itime_in[2], int n_lats, int n_lons, int model_number, string filename) {
    // Run the field-line tracer and output the results. Cute!

    FILE * file;

    double lmin = 10;
    double lmax = 90;

    double lat_spacing = (lmax - lmin)/(1.0*n_lats);
    double lon_spacing = 360./(1.0*n_lons);

    int n_rays = n_lats * n_lons;
    double grid[n_rays][TRACER_MAX][3];
    double lens[n_rays]; 

    double lats[n_rays];
    double lons[n_rays];

    double x_geo[3], x_geocar[3], x_sm[3];

    double tsyg_params[10] = {0};
    double VG[3];
    double stepsize = 0.01;

    // Setup for IGRF:    
    load_TS05_params(itime_in, tsyg_params, VG);
    init_igrf(itime_in);

    // set DST:
    tsyg_params[1] = -20;


    // Trace it
    int i = 0;
    for (double lat = lmin; lat < lmax; lat+= lat_spacing) {
        for (double lon = 0; lon < 360; lon += lon_spacing) {
            cout << "tracing " << lat << ", " << lon << "\n";

            init_igrf(itime_in);
            lats[i] = lat;
            lons[i] = lon;

            x_geo = {1.2, lat, lon};
            x_geocar = {x_geo[0], x_geo[1], x_geo[2]};
            degcar(x_geocar);
            mag_to_sm_d_(itime_in, x_geocar, x_sm);

            lens[i] = trace_fieldline(itime_in, x_sm, grid[i], NULL, stepsize, model_number, tsyg_params);
            i++;
        }
    }


    // Save it
    cout << "saving to file " << filename << "\n";
    file = fopen(filename.c_str(), "w");

    if (file != NULL) {
        int k = 0;
        for (int i=0; i < n_rays; i++) {
            cout << "Saving " << lats[k] << ", " << lons[k] <<"\n";
            for (int j=0; j < lens[i]; j++) {
                // cout << "(" << lats[k] << ", " << lons[k] << "): " << j << " ";
                // print_array(grid[i][j], 3);

                fprintf(file, "%g %g %i %g %g %g\n",lats[k], lons[k], j, grid[i][j][0], grid[i][j][1], grid[i][j][2]);

            }
            k++;
            // cout << "\n";
        }
    } else {
        cout << "Could not open file " << filename.c_str() << "\n";
    }

}



bool coarse_mask(rayT cur_frames[8], rayT prev_frames[8], EA_segment EA) {
    // Coarse masking detection. Returns false if all ray points are
    // on the same side of the EA_arr plane. True otherwise.

    Vector3d EApos = EA.ea_pos;
    Vector3d n = EA.ea_norm;
    Vector3d l1, l2;
    double s1, s2, sr1, sr2, sl1, sl2, rad1, rad2;

    double EAr[3];
    double r1[3], r2[3];

    double side = 0;
    double rside= 0;
    double lside = 0;

    bool sideflag = false;
    bool radflag  = false;
    bool lonflag  = false;

    carsph(EA.ea_pos.data(), EAr);
    
    double lon_width = ((EA.radius)/EAr[0]);
    // cout << "width: " << R2D*lon_width << "\n";

    // cout << EAr[0] << ", " << EAr[1]*R2D << ", " << EAr[2]*R2D << "\n";

    for (int rr=0; rr < 8; rr++) {
        // l1 = Map<VectorXd>(cur_rays[rr]->pos[t-1].data(), 3,  1);
        // l2 = Map<VectorXd>(cur_rays[rr]->pos[t  ].data(), 3,  1);

        l1 = prev_frames[rr].pos;
        l2 = cur_frames[rr].pos;
        carsph(l1.data(), r1);
        carsph(l2.data(), r2);

        // Check if any are on opposite sides of the plane:
        if (rr==0) {
            side = signof((l1 - EApos).dot(n));
        }

        s1 = signof((l1 - EApos).dot(n));
        s2 = signof((l2 - EApos).dot(n));

        if ( (s1 != side) || (s2 != side) ) {
            sideflag = true;
            // return true;
        }

        // Check if radii are all too low or too high:
        if (rr==0) {
            rside= ( (r1[0] < EAr[0] - L_MARGIN) ? - 1 :
                     (r1[0] > EAr[0] + L_MARGIN) ?   1 : 0);
        }
        sr1 = ( (r1[0] < EAr[0] - L_MARGIN) ? - 1 :
                (r1[0] > EAr[0] + L_MARGIN) ?   1 : 0);

        sr2 = ( (r2[0] < EAr[0] - L_MARGIN) ? - 1 :
                (r2[0] > EAr[0] + L_MARGIN) ?   1 : 0);

        if ( (sr1 != rside) || (sr2 != rside) ) {
            radflag = true;
        }

        // cout << "lon1: " << R2D*r1[2] << " EAlon: " << R2D*EAr[2] << "\n";


        // Check if longitudes are all on the same side of the EA:

        if (rr==0) {
            lside = longitude_interval(r1[2], EAr[2], lon_width);
        }
        sl1 = longitude_interval(r1[2], EAr[2], lon_width);
        sl2 = longitude_interval(r2[2], EAr[2], lon_width);

        if ( (sl1 != lside) || (sl2 != lside) ) {
            lonflag = true;
        }
    }

    // All radii are between the limits!
    if (rside == 0) {
        radflag = true;
    }
    // All longitudes in sweet zone!
    if (lside == 0) {
        // cout << "derp ";
        lonflag = true;
    }
    return sideflag && radflag && lonflag; 

}


double longitude_interval(double ra, double r0, double width) {
    // Pretty sure this won't wrap around happily, but it works for now.
    // (10.28.16)

    // cout << "ra: " << ra*R2D << " r0: " << r0*R2D << " diff: " << R2D*(ra - r0) << "\n";

    // const double width = 5;

    if ( (ra - r0) < - width) {
        return -1;
    } else if ( (ra - r0) > width) {
        return 1;
    } else {
        return 0;
    }
}


bool crosses_EA(Vector3d p1, Vector3d p2, EA_segment EA) {
    // Returns true if the vector between points p1 and p2 passes
    // through the EA segment EA (specified by point p0 and normal n)
    // see: http://paulbourke.net/geometry/pointlineplane/

    double u, rr;

    Vector3d p0 = EA.ea_pos;
    Vector3d n  = EA.ea_norm;
    Vector3d p;                 // intersection point

    u = ( (p0 - p1).dot(n) ) / ( (p2 - p1).dot(n) );
    p = p1 + u*(p2 - p1);
    rr = (p - p0).norm();

    if (abs(u) <= 1) {
        return (rr <= EA.radius);
    } else {
        return false;
    }

}


void dump_EA_array(vector<EA_segment> EA_array, string filename) {
    // Write the EA array to a file, so I can plot it.

    FILE * file;

        // Save it
    cout << "saving EA segments to file " << filename << "\n";
    file = fopen(filename.c_str(), "w");
    if (file != NULL) {
        
        for (int i=0; i < EA_array.size(); i++) {
            fprintf(file, "%g %g %g %g %g %g %g\n",
                EA_array[i].ea_pos[0], EA_array[i].ea_pos[1],  EA_array[i].ea_pos[2],
                EA_array[i].ea_norm[0],EA_array[i].ea_norm[1], EA_array[i].ea_norm[2],
                EA_array[i].radius);
        }
    fclose(file);
    } else {
        cout << "Could not open file " << filename.c_str() << "\n";
    }
}


void add_rayT(rayT* rayA, rayT* rayB) {

    rayA->inp_pwr += rayB->inp_pwr;       
    rayA->damping += rayB->damping;      
    rayA->n       += rayB->n;         
    rayA->B0      += rayB->B0;        
    rayA->stixP   += rayB->stixP;         
    rayA->stixR   += rayB->stixR;         
    rayA->stixL   += rayB->stixL;         
        
    // rayA->stixS   += rayB->stixS;         
    // rayA->stixD   += rayB->stixD;

    rayA->ds      += rayB->ds;

    rayA->num_rays += 1;         

    for (int i=0; i < rayA->Ns.size(); i++) {
        rayA->Ns[i]  += rayB->Ns[i];
        rayA->nus[i] += rayB->nus[i];        
    }
}


void interp_rayF(rayF* rayfile, rayT* frame, double t_target) {
// float interpPt(float *xI, float *yI, int n, float t_target )
                //  time axis   data vector length   t_targ       
    int i, iHigh, iLow, iMid;
    double yO;

    double M;
    vector <double> xI = rayfile->time;
    int n = rayfile->time.size();


    // Check that t_target is within bounds
    if( (t_target < xI[0]) || t_target  > xI[n-1] ) {
        printf("\nPoint is out of bounds! %g, {%g, %g}\a\n",t_target , xI[0], xI[n-1]);
        return;
    }
      
    // Do a binary search for the correct index 
    iHigh = n-1;
    iLow = 0;  
    while(1) {
        iMid = (iHigh+iLow)/2;
        if( (t_target  >= xI[iLow]) && (t_target  < xI[iMid]) ) {
            iHigh = iMid;
        } else {
            iLow = iMid;
        }
        if(t_target ==xI[n-1]){printf("\nin interpPt\n"); return;}
        if(iHigh==iLow) {
            printf("\nexiting from binary search in 1st condtion\n");
            break;
        }
        if( (t_target >=xI[iMid]) && (t_target <xI[iMid+1]) ) break;
        if( (t_target >=xI[iMid-1]) && (t_target <xI[iMid]) ) {
           iMid--;
           break;
            }
        }

    M = ( t_target -xI[iMid] )/( xI[iMid+1]-xI[iMid] );
    // Now, let's interpolate the output values:
    // Vector-valued
    for (int k = 0; k < 3; k++) {
        frame->pos[k] = ( rayfile->pos[iMid+1][k]-rayfile->pos[iMid][k] )*M + rayfile->pos[iMid][k];
        frame->n[k] =   ( rayfile->n[iMid+1][k]  -rayfile->n[iMid][k]   )*M + rayfile->n[iMid][k];
        frame->B0[k] =  ( rayfile->B0[iMid+1][k] -rayfile->B0[iMid][k]  )*M + rayfile->B0[iMid][k];
    }
    // Scalar-valued
    frame->damping = ( rayfile->damping[iMid+1]-rayfile->damping[iMid] )*M + rayfile->damping[iMid];
    frame->stixP = ( rayfile->stixP[iMid+1]-rayfile->stixP[iMid] )*M + rayfile->stixP[iMid];
    frame->stixR = ( rayfile->stixR[iMid+1]-rayfile->stixR[iMid] )*M + rayfile->stixR[iMid];
    frame->stixL = ( rayfile->stixL[iMid+1]-rayfile->stixL[iMid] )*M + rayfile->stixL[iMid];
    
    // Stuff that doesn't need interpolation:
    frame->w    = rayfile->w;
    frame->in_lat = rayfile->in_lat;
    frame->in_lon = rayfile->in_lon;
    frame->time = t_target;
    frame->inp_pwr = rayfile->inp_pwr;
}


// vector <vector <int> > find_adjacent_rays(map <int, vector<double> > start_locs) {
//     // Given a set of input starting coordinates, determine
//     // a list of adjacent rays to iterate over.
//     // (this is kind of overkill)

//     // map <int, int> closest_up;
//     map <int, int> closest_down;
//     map <int, int> closest_right;

//     double lat1, lat2, lon1, lon2;
//     double dist_down, dist_right, dist;
//     int cur_ind, check_ind;
//     double lat3, lon3, lat4, lon4;
//     double d1, d2;
    
//     vector <vector <int> > adjacency_list;


//     cout << "Finding adjacent rays...\n";
//     for(map<int, vector<double> >::iterator iter = start_locs.begin(); iter != start_locs.end(); ++iter){
//         // cout << "searching around ";  print_vector(iter->second);
//         cur_ind = iter->first;

//         lat1 = iter->second[1];
//         lon1 = iter->second[2];
//         dist_down = 1e12;
//         dist_right= 1e12;

//         for(map<int, vector<double> >::iterator other = start_locs.begin(); other != start_locs.end(); ++other){
//             if (iter->first != other->first) {
//                 lat2 = other->second[1];
//                 lon2 = other->second[2];
//                 check_ind = other->first;

//                 dist = haversine_distance(lat1, lon1, lat2, lon2);

//                 // find closest entry below current latitude
//                 if ( lat1 - lat2 > 0.001) {
//                     // Distance in latitude:
//                     // printf("dist from (%2.1f, %2.1f) to (%2.1f, %2.1f): %g\n",lat1, lon1, lat2, lon2, dist);
//                     if (dist < dist_down) {
//                         closest_down[iter->first] = other->first;
//                         dist_down = dist;
//                         // cout << "dist_down = " << dist_down << " ";
//                         // cout << "iter: " << iter->first << " other: " << other->first << "\n";
//                     }
//                 }
//                 // find closest entry to right of current longitude
//                 if (lon2 - lon1 > 0.001) {
//                     if (dist < dist_right) {
//                         closest_right[iter->first] = other->first;
//                         dist_right = dist;
//                     }
//                 }
//             }   // not checking against itself
//         }   // Loop over each entry
//     }



//     // Next, select entries which have both a ray down and a ray right:
//     for(map<int, vector<double> >::iterator iter = start_locs.begin(); iter != start_locs.end(); ++iter){
//         vector <int> cur_inds;

//         // cout << iter->first << " ";
//         if ( (closest_down.count(iter->first) != 0) && (closest_right.count(iter->first) != 0) ) {
//             lat1 = iter->second[1];
//             lon1 = iter->second[2];
//             lat2 = start_locs.at(closest_down[iter->first])[1];
//             lon2 = start_locs.at(closest_down[iter->first])[2];
//             lat3 = start_locs.at(closest_right[iter->first])[1];
//             lon3 = start_locs.at(closest_right[iter->first])[2];
        
//             // Select the fourth value -- is it right and down, or down and right?
//             d1 = 1e12;
//             d2 = 1e12;

//             if (closest_down.count(closest_right[iter->first]) != 0) {
//                 lat4 = start_locs[closest_down[closest_right[iter->first]]][1];
//                 lon4 = start_locs[closest_down[closest_right[iter->first]]][2];
//                 d1 = haversine_distance(lat1, lon1, lat4, lon4);
//             }
//                 // cout << "down and right: " << closest_down[closest_right[iter->first]] << "\n";
//             if (closest_right.count(closest_down[iter->first]) != 0) {
//                 lat4 = start_locs[closest_right[closest_down[iter->first]]][1];
//                 lon4 = start_locs[closest_right[closest_down[iter->first]]][2];
//                 d2 = haversine_distance(lat1, lon1, lat4, lon4);

//                 // cout << "right and down: " << closest_right[closest_down[iter->first]] << "\n";
//             }

//             if ( (d1 != 1e12) || (d2 != 1e12) ) {
//                 cur_inds.push_back(iter->first);
//                 cur_inds.push_back(closest_right[iter->first]);
//                 cur_inds.push_back(closest_down[iter->first]);
//                 if (d1 <= d2) {
//                     cur_inds.push_back(closest_down[closest_right[iter->first]]);
//                 } else {
//                     cur_inds.push_back(closest_right[closest_down[iter->first]]);
//                 }
//                 adjacency_list.push_back(cur_inds);
//             }
//             // printf("ray at (%2.1f, %2.1f)  -> (%2.1f, %2.1f), (%2.1f, %2.1f)\n",lat1, lon1, lat2, lon2, lat3, lon3);
//         }
//     }

//     // cout << "Found " << adjacency_list.size() << " sets of adjacent guide rays\n";
//     // for (int i=0; i < adjacency_list.size(); i++ ){

//     //     cout << adjacency_list[i][0] << ", ";
//     //     cout << adjacency_list[i][1] << ", ";
//     //     cout << adjacency_list[i][2] << ", ";
//     //     cout << adjacency_list[i][3] << "\n";
//     // }

//     // // Print out some results:
//     // cout << "closest down: \n";
//     // for(map<int, int>::iterator iter = closest_down.begin(); iter != closest_down.end(); ++iter){
//     //     lat1 = start_locs.at(iter->first)[1];
//     //     lon1 = start_locs.at(iter->first)[2];
//     //     lat2 = start_locs.at(iter->second)[1];
//     //     lon2 = start_locs.at(iter->second)[2];
//     //     cout << iter->first << " " << iter->second << " ";
//     //     printf("ray at (%2.1f, %2.1f)  -> (%2.1f, %2.1f)\n",lat1, lon1, lat2, lon2);

//     // }
//     // cout << "closest right: \n";
//     // for(map<int, int>::iterator iter = closest_right.begin(); iter != closest_right.end(); ++iter){
//     //     lat1 = start_locs.at(iter->first)[1];
//     //     lon1 = start_locs.at(iter->first)[2];
//     //     lat2 = start_locs.at(iter->second)[1];
//     //     lon2 = start_locs.at(iter->second)[2];
//     //     cout << iter->first << " " << iter->second << " ";
//     //     printf("ray at (%2.1f, %2.1f)  -> (%2.1f, %2.1f)\n",lat1, lon1, lat2, lon2);

//     // }

// return adjacency_list;

// }


vector< vector<double> > find_adjacent_rays(vector< vector<double> > available_rays) {
    // Returns a list of four (lat, lon) pairs, corresponding to adjacent rays
    vector< vector<double> > adjacent_rays;
    vector<double>::iterator it;
    vector<double> uLats, uLons;
    double row[8];

    for (vector< vector<double> >::iterator itt=available_rays.begin(); itt!=available_rays.end(); ++itt) {
        // print_vector(*itt);
        uLats.push_back((*itt)[1]);
        uLons.push_back((*itt)[2]);
    }

    // Sorted list of unique latitudes
    // copy(start_lats.begin(), start_lats.end(), back_inserter(uLats));
    sort(uLats.begin(), uLats.end());
    it = unique(uLats.begin(), uLats.end());
    uLats.resize(distance(uLats.begin(), it));
    cout << "ray lats: ";
    print_vector(uLats);
    // Sorted list of unique longitudes
    // copy(start_lons.begin(), start_lons.end(), back_inserter(uLons));
    sort(uLons.begin(), uLons.end());
    it = unique(uLons.begin(), uLons.end());
    uLons.resize(distance(uLons.begin(), it));

    cout << "ray lons: ";
    print_vector(uLons);

    for (int la = 0; la < uLats.size()-1; ++la) {
        for (int lo = 0; lo < uLons.size()-1; ++lo) {
              
            row = {uLats[la],   uLons[lo],   
                   uLats[la+1], uLons[lo],   
                   uLats[la],   uLons[lo+1], 
                   uLats[la+1], uLons[lo+1]};

            adjacent_rays.push_back(vector<double>(row, row+8));
        }
    }

    return adjacent_rays;
}


vector< vector<double> > find_adjacent_rays_2d(vector< vector<double> > available_rays) {
    // Returns a list of latitude pairs, corresponding to adjacent rays
    vector< vector<double> > adjacent_rays;
    vector<double>::iterator it;
    vector<double> uLats, uLons;
    double row[8];

    for (vector< vector<double> >::iterator itt=available_rays.begin(); itt!=available_rays.end(); ++itt) {
        // print_vector(*itt);
        uLats.push_back((*itt)[1]);
        uLons.push_back((*itt)[2]);
    }

    // Sorted list of unique latitudes
    // copy(start_lats.begin(), start_lats.end(), back_inserter(uLats));
    sort(uLats.begin(), uLats.end());
    it = unique(uLats.begin(), uLats.end());
    uLats.resize(distance(uLats.begin(), it));
    cout << "ray lats: ";
    print_vector(uLats);
    // Sorted list of unique longitudes
    // copy(start_lons.begin(), start_lons.end(), back_inserter(uLons));
    sort(uLons.begin(), uLons.end());
    it = unique(uLons.begin(), uLons.end());
    uLons.resize(distance(uLons.begin(), it));

    cout << "ray lons: ";
    print_vector(uLons);

    for (int la = 0; la < uLats.size()-1; ++la) {
        for (int lo = 0; lo < uLons.size(); ++lo) {
              
            row = {uLats[la],   uLons[lo],   
                   uLats[la+1], uLons[lo],
                   uLats[la],   uLons[lo],   
                   uLats[la+1], uLons[lo]};

            adjacent_rays.push_back(vector<double>(row, row+8));
        }
    }

    return adjacent_rays;
}

















cellT new_cell(rayT ray) {
    cellT cell;

    cell.pos = ray.pos;
    cell.t   = ray.time;
    cell.f   = ray.w/(2*PI);
    // cell.pwr = ray.inp_pwr*ray.damping/FREQ_STEP_SIZE/TIME_STEP/ray.ds;
    // cell.pwr = ray.inp_pwr*ray.damping;
    // cell.pwr = ray.inp_pwr;

    Vector3d kvec = ray.n*ray.w/C;
    double k      = kvec.norm();
    Vector3d Bhat = ray.B0.array()/ray.B0.norm();
    double kz     = -1.0*kvec.dot(Bhat); //k.array()*Bhat.array();
    double kx     = (kvec + kz*Bhat).norm();

    // angle between parallel and perpendicular K
    cell.psi = R2D*atan2(-kx, kz);

    // cout << "psi: " << cell.psi << "\n";
    cell.mu = ray.n.norm();
    cell.stixP    = ray.stixP;
    cell.stixR    = ray.stixR;
    cell.stixL    = ray.stixL;
    cell.num_rays = 1;

    return cell;
}


void add_cell(cellT* cell1, cellT* cell2) {
    
    cell1->pwr       += cell2->pwr;
    cell1->psi       += cell2->psi;
    cell1->mu        += cell2->mu;
    cell1->stixP     += cell2->stixP;
    cell1->stixR     += cell2->stixR;
    cell1->stixL     += cell2->stixL;
    cell1->num_rays  += 1; 

    for (int i=0; i< cell1->pwr_vec.size(); i++) {
        cell1->pwr_vec[i] += cell2->pwr_vec[i];
    }
}



void calc_resonance(map<pair<int,int>, cellT> db, EA_segment EA, 
    double da_N[NUM_E][NUM_TIMES], double da_S[NUM_E][NUM_TIMES], int lon_index) {

    int i, j=0, kk, mres, noutN, noutS, ei, ti, e_toti;
    cellT *next;
    // char *prefN="pN", *prefS="pS", suff[64];
    double lat, L, t, f, pwr, psi, mu, stixP, stixR, stixL, latk;
    double Bxw, Byw, Bzw, Exw, Eyw, Ezw, stixD, stixS, stixA;
    double stixB, stixX, n_x, n_z, k, kx, kz, rho1, rho2, Byw_sq;
    double flt_const_N, flt_const_S, flt_time, eta_dot;

    // double flt_const_N[EA_SPLIT], flt_const_S[EA_SPLIT], flt_time, eta_dot;
    double wh, dwh_ds, gamma, alpha1, alpha2, beta, v_para, v_perp;
    double spsi, cpsi, spsi_sq, cpsi_sq, mu_sq, w, R1, R2, w1, w2;
    double alpha_lc, alpha_eq, epsm, slat, clat, slat_term, ds;
    double t1, t2, t3, direction, v_para_res, v_tot, v_tot_res;
    double salph, calph, wtau_sq, Y, dv_para_ds, AA, BB, T1;
    double Farg, Farg0, Fs, Fc, Fs0, Fc0, dFs_sq, dFc_sq, dalpha;
    double alpha_eq_p, dalpha_eq, E_res, e_starti, e_endi;
    double v_tot_arr[NUM_E], E_tot_arr[NUM_E], v_para_star, v_para_star_sq;
    float *arr_N, *arr_S;
    time_t start, end;
    int timei;
    cellT cell;
  
    L = EA.Lsh;
    epsm = (1/L)*(R_E+H_IONO)/R_E;
    // alpha_eq = asin(sqrt( pow(epsm,3)/sqrt(1+3*(1-epsm)) ));
    alpha_eq = EA.alpha_eq;
    //initialize the velocity and energy arrays
    for(i=0; i<NUM_E; i++) {
        E_tot_arr[i] = pow(10, (E_EXP_BOT+ DE_EXP*i) ); // energy in eV
        v_tot_arr[i] = C*sqrt(1 - pow( (E_EL/(E_EL+E_tot_arr[i])) ,2) );
    }

    // print_array(v_tot_arr, 10);
    
    // cout << "v_tot_arr[0]: " << v_tot_arr[0] << " v_tot_arr[end]: " << v_tot_arr[NUM_E-1] << "\n";
    // cout << "E_tot_arr[0]: " << E_tot_arr[0] << " E_tot_arr[end]: " << E_tot_arr[NUM_E-1] << "\n";
    lat = EA.lat;
    slat = sin( lat*D2R );
    clat = cos( lat*D2R );
    slat_term = sqrt(1+3*slat*slat);
    // wh = 2*PI*880000/pow(L,3)*slat_term/pow(clat,6);
    // dwh_ds = 3*wh/(L*R_E)*slat/slat_term*
      // (1/(slat_term*slat_term) + 2/(clat*clat));
    wh = EA.wh;
    dwh_ds = EA.dwh_ds;


    // getFltConst(L,lat,alpha_eq,&(flt_const_N),&(flt_const_S));
    flt_const_N = EA.ftc_n;
    flt_const_S = EA.ftc_s;

    alpha_lc = EA.alpha_lc; 
    // alpha_lc = asin(sqrt( slat_term/pow(clat,6) )*sin(alpha_eq));
    salph = sin(alpha_lc);
    calph = cos(alpha_lc);
    ds = EA.ds;
    dv_para_ds = EA.dv_para_ds;

    printf("------------- EA at lat: %2.2f ---------------\n", lat);
    
    // Go through the cells in each latitude
    for (map<pair<int,int>,cellT>::iterator iter = db.begin(); iter != db.end(); iter++) {
        cell = iter->second;

        t = cell.t + TIME_STEP/2;           // We want the time and freq to be in the 
        f = cell.f + FREQ_STEP_SIZE/2;      // center of the cell, so add DT/2 or DF/2



        // pwr = cell.pwr*FREQ_STEP_SIZE/cell.num_rays;
        // pwr = cell.pwr/cell.num_rays/TIME_STEP/EA.ds;
        // pwr = sqrt(cell.pwr)/EA.ds;  
        // pwr = sqrt(cell.pwr);

        // cout << "lon_index = " << lon_index << endl;
        // cout << "size of power vector: " << cell.pwr_vec.size() << endl;
        // print_vector(cell.pwr_vec);

        pwr = cell.pwr_vec[lon_index]/cell.num_rays;
        // pwr = cell.pwr/cell.num_rays;

        psi = D2R*cell.psi/cell.num_rays;

        printf("t: %g f: %g pwr: %g \n",t,f,pwr);

        if(pwr > WAVE_PWR_THRESH) {


            mu = cell.mu/cell.num_rays;
            stixP = cell.stixP/cell.num_rays;
            stixR = cell.stixR/cell.num_rays;
            stixL = cell.stixL/cell.num_rays;

            spsi = sin(psi);
            cpsi = cos(psi);
            spsi_sq = pow(spsi,2);
            cpsi_sq = pow(cpsi,2);
            n_x = mu*fabs(spsi);
            n_z = mu*cpsi;
            mu_sq = mu*mu;
            w = 2.0*PI*f;
            k = w*mu/C;
            kx = w*n_x/C;
            kz = w*n_z/C;
            Y = wh / w ;

            // Stix parameters
            stixS = ( stixR + stixL ) /2.0;
            stixD = ( stixR - stixL ) /2.0;
            stixA = stixS + (stixP-stixS)*cpsi_sq;
            stixB = stixP*stixS+stixR*stixL+(stixP*stixS-stixR*stixL)*cpsi_sq;
            stixX = stixP/(stixP- mu_sq*spsi_sq);

            // Polarization ratios
            rho1=((mu_sq-stixS)*mu_sq*spsi*cpsi)/(stixD*(mu_sq*spsi_sq-stixP));
            rho2 = (mu_sq - stixS) / stixD ;

            // (bortnik 2.28)
            Byw_sq =  2.0*MU0/C*pwr*stixX*stixX*rho2*rho2*mu*fabs(cpsi)/
               sqrt(  pow((tan(psi)-rho1*rho2*stixX),2) + 
               pow( (1+rho2*rho2*stixX), 2 ) );


            // RMS wave components
            Byw = sqrt(Byw_sq);
            Exw = fabs(C*Byw * (stixP - n_x*n_x)/(stixP*n_z)); 
            Eyw = fabs(Exw * stixD/(stixS-mu_sq));
            Ezw = fabs(Exw *n_x*n_z / (n_x*n_x - stixP));
            Bxw = fabs(Exw *stixD*n_z /C/ (stixS - mu_sq));
            Bzw = fabs((Exw *stixD *n_x) /(C*(stixX - mu_sq)));

            // Oblique integration quantities
            R1 = (Exw + Eyw)/(Bxw+Byw);
            R2 = (Exw - Eyw)/(Bxw-Byw);
            w1 = Q_EL/(2*M_EL)*(Bxw+Byw);
            w2 = Q_EL/(2*M_EL)*(Bxw-Byw);
            alpha1 = w2/w1;

            if (DEBUG) {
                printf("Byw: %g Exw: %g Eyw: %g Ezw: %g Bxw: %g Bzw: %g\n",
                        Byw,    Exw,    Eyw,    Ezw,    Bxw,    Bzw);
                printf("t: %g, f: %g, pwr: %g, psi: %g, Num_rays: %g\n",t,f,pwr, R2D*psi, cell.num_rays);
                printf("wh: %g dwh_ds: %g alpha_lc: %g alpha_eq: %g ds: %g dv_para_ds: %g\n",
                        wh,    dwh_ds,    alpha_lc,    alpha_eq,    ds,    dv_para_ds);
                printf("R1: %g R2: %g w1: %g w2: %g alpha1: %g\n",
                        R1,    R2,    w1,    w2,    alpha1);   
            }


            //begin MRES loop here
            for(mres=-SCATTERING_RES_MODES; mres <= SCATTERING_RES_MODES; mres++) {
                // get parallel resonance velocity
                t1 = w*w*kz*kz;
                t2 = pow((mres*wh),2)-w*w;
                t3 = kz*kz + pow((mres*wh),2)/(pow(C*cos(alpha_lc),2));

                if(mres==0) {
                    direction = -kz/fabs(kz);
                } else {
                    direction = kz/fabs(kz) * mres/fabs(mres) ;
                }

                v_para_res = ( direction*sqrt(t1 + t2*t3) - w*kz ) / t3;
                v_tot_res = v_para_res / cos(alpha_lc); 
                E_res = E_EL*( 1.0/sqrt( 1.0-(v_tot_res*v_tot_res/(C*C)) ) -1.0 );

                // if(DEBUG) {printf("t1: %g t2: %g t3: %g v_para_res: %g v_tot_res: %g E_res: %g\n",
                //                    t1,    t2,    t3,    v_para_res,    v_tot_res,    E_res);}
                
                // get starting and ending indices, +-20% energy band
                e_starti = floor((log10(E_res) - E_EXP_BOT - 0.3)/(DE_EXP));
                e_endi   =  ceil((log10(E_res) - E_EXP_BOT + 0.3)/(DE_EXP));

                if(e_endi>NUM_E) e_endi=NUM_E;
                if(e_starti>NUM_E) e_starti=NUM_E;
                if(e_endi<0) e_endi=0;
                if(e_starti<0) e_starti=0;
                

                // begin V_TOT loop here
                for(e_toti=e_starti; e_toti < e_endi; e_toti++) {

                    v_tot = direction*v_tot_arr[e_toti];
                    v_para = v_tot * calph;
                    v_perp = fabs(v_tot * salph);

                    gamma = 1.0 / sqrt(1 - pow((v_tot/C),2)); 
                    alpha2 = Q_EL*Ezw /(M_EL*gamma*w1*v_perp);
                    beta = kx*v_perp / wh ;
                    wtau_sq = pow((-1),(mres-1)) * w1/gamma * 
                    ( jn( (mres-1), beta ) - 
                      alpha1*jn( (mres+1) , beta ) +
                      gamma*alpha2*jn( mres , beta ) ); 
                    T1 = -wtau_sq*(1+ ( (calph*calph) / (mres*Y-1) ) );

                    // Now - start analytical evaluation!!!
                  
                    if( fabs(lat)< 1e-3) {
                        // Near the equator we can use a simplified expression:
                        eta_dot = mres*wh/gamma - w - kz*v_para;

                        if(fabs(eta_dot)<10) {
                            // Bortnik A.31
                            dalpha_eq = fabs(T1/v_para)*ds/sqrt(2); 
                        } else {
                            // Bortnik A.30
                            dalpha_eq = fabs(T1/eta_dot)*sqrt(1-cos(ds*eta_dot/v_para)); 
                        }

                    } else {  
                        
                        v_para_star = v_para - dv_para_ds*ds/2.0;
                        v_para_star_sq = v_para_star * v_para_star;

                        // Bortnik A.18 -- part A1
                        AA = (mres/(2.0*v_para_star*gamma))*dwh_ds* 
                             (1 + ds/(2.0*v_para_star)*dv_para_ds) - 
                              mres/(2.0*v_para_star_sq*gamma)*wh*dv_para_ds + 
                              w/(2.0*v_para_star_sq)*dv_para_ds ;

                        // // Bortnik A.18 -- part A0   -- THIS DOES NOT MATCH THE THESIS
                        // BB = mres/(gamma*v_para_star)*wh - 
                        //      mres/(gamma*v_para_star)*dwh_ds*(ds/2.0) -
                        //      w/v_para_star - kz;

                        // Bortnik A.18 -- part A0
                        BB =   mres*wh/(gamma*v_para_star)
                             - mres/(gamma*v_para_star)*dwh_ds*(ds/2.0) * (w/v_para_star)*kz;


                        // Evaluate Bortnik A.26 -- integration performed thru Fresnel functions
                        Farg = (BB + 2*AA*ds) / sqrt(2*PI*fabs(AA));
                        Farg0 = BB / sqrt(2*PI*fabs(AA));  
                        
                        Fresnel(Farg, &Fs, &Fc);
                        Fresnel(Farg0, &Fs0, &Fc0);
                        
                        dFs_sq = pow((Fs - Fs0),2);
                        dFc_sq = pow((Fc - Fc0),2);
                        
                        dalpha = sqrt(PI/4/fabs(AA))*fabs(T1/v_para)*sqrt(dFs_sq+dFc_sq);
                        
                        // Map the local change in pitch angle to the equivalent
                        // pitch angle at the equator:  (still using dipole model here)
                        // alpha_eq_p = asin( sin(alpha_lc+dalpha)*pow(clat,3) / 
                        //            sqrt(slat_term) );
                        alpha_eq_p = asin( sin(alpha_lc + dalpha)*EA.Bo_ratio);

                        dalpha_eq = alpha_eq_p - alpha_eq;

                        
                        if (isnan(dalpha_eq)) {
                            cout << "NaN: ";

                            // cout << "direction: " << direction << " v_to_arr[e_toti]: " << v_tot_arr[e_toti] << " "; 
                            // printf("w1: %g gamma: %g alpha1: %g beta: %g alpha2: %g v_perp: %g v_tot: %g salph: %g e_toti: %g\n",
                            //         w1,    gamma,    alpha1,    beta,    alpha2,    v_perp,    v_tot,    salph,    e_toti);


                            // printf("wtau_sq: %g calph: %g mres: %g Y: %g\n",
                            //         wtau_sq,    calph,    mres,    Y);
                            // printf("AA: %g T1: %g v_para: %g dFs_sq: %g dFc_sq: %g\n",
                            //         AA,     T1,    v_para,   dFs_sq,    dFc_sq);
                            // // printf("alpha_eq_p: %g alpha_lc: %g dalpha: %g Bo_ratio: %g\n",
                            //         alpha_eq_p,    alpha_lc,    dalpha,  EA.Bo_ratio);
                            // printf("flt_time: %g dalpha_eq: %g alpha_eq_p: %g alpha_eq: %g\n",
                            //         flt_time,    dalpha_eq,    alpha_eq_p,    alpha_eq);
                            
                            break;
                        }
                    }

                    if(direction>0) {
                        flt_time = fabs(flt_const_N/v_para);
                    } else {
                        flt_time = fabs(flt_const_S/v_para);
                    }
                     
                    // if (DEBUG) {printf("flt_time: %g dalpha_eq: %g alpha_eq_p: %g alpha_eq: %g\n",
                    //                     flt_time,    dalpha_eq,    alpha_eq_p,    alpha_eq);}

                    // Get time index into output array
                    timei = round((t + flt_time)/TIME_STEP);

                    if (timei < NUM_TIMES) {
                        // Save it!
                        if (direction > 0) {
                            da_N[e_toti][timei] += dalpha_eq*dalpha_eq;
                        } else {
                            da_S[e_toti][timei] += dalpha_eq*dalpha_eq;
                        }
                    } else {
                        // Do we want to track total scattering after TMAX?
                    }
                } // v_para
            } // mres loop
        } // if pwr > 0.0
    } // entries in db
}





double polygon_frame_area(rayT frame[8]) {
    // Calculates the area, in meters, enclosed by the set of guide rays.

    Vector3d cp(0,0,0);
    int inds[4] = {0,1,2,3};
    int n=4;
    
    double max_area = 0;
    double area;

    do {
        area = 0;
        Vector3d cp(0,0,0);

        for (int i=0; i<n; ++i) {
            Vector3d v1 = frame[inds[i]].pos;
            Vector3d v2 = frame[inds[(i+1)%n]].pos;
            cp += v1.cross(v2);
        }

        area = pow(R_E, 2)*cp.norm()/2.;

        if (area > max_area) { max_area = area;}
    } while ( next_permutation(inds,inds + 2) );

    return max_area;
}

