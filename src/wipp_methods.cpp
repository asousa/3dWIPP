#include <wipp.h>

using namespace std;
using namespace Eigen;


double input_power_scaling(double* flash_loc, double* ray_loc, double mag_lat, double w, double i0) {
    // Returns ray power at the top of the ionosphere
    // per unit in space and frequency.

    double theta;        // angle between two vectors
    double gc_distance; // great circle distance between two points
    double dist_tot;
    double xi;

    double S, S_vert;
    double attn_factor;
    double w_sq, f;

    Vector3d v1;
    Vector3d v2;


    f = w/(2.0*PI);

    // cout << "flash loc: " << flash_loc[0] << " " << flash_loc[1] << " " << flash_loc[2] << "\n";
    // cout << "flash loc: " << ray_loc[0] << " " << ray_loc[1] << " " << ray_loc[2] << "\n";
   
   
    v1 = Map<VectorXd>(flash_loc,3,1);
    v2 = Map<VectorXd>(ray_loc,  3,1);

    theta = acos(v1.dot(v2)/(v1.norm()*v2.norm()));
    
    // cout << "theta: " << R2D*theta << " deg\n";
    // Arc length (~great-circle distance) between vI0ectors
    gc_distance = (R_E)*theta;
    // cout << " gc_dist: " << gc_distance*1e-3 << "\n";
    // total distance up to ionosphere:
    dist_tot = hypot(gc_distance, H_IONO);
    xi = atan2(gc_distance, H_IONO);  // Incident angle

    // cout << "xi: " << R2D*xi << " deg\n";
    w_sq =  pow( w , 2 );
    S = ( (1/Z0) * pow( (H_E*i0*2E-7*(sin(xi)/dist_tot)*w*(P_A-P_B)) , 2 ) 
                   /  (  (w_sq+pow(P_A,2))*(w_sq+pow(P_B,2))  )      ) ;
    S_vert = S * cos(xi) ;  // factor for vert prop.


    // Ionosphere absorption model
    attn_factor = pow(10,-(ionoAbsorp(mag_lat,f)/10)  );
    S_vert = S_vert * attn_factor;

    printf("i0: %2.3f, dist_tot: %2.3f, xi: %2.3f, S_vert; %e\n",i0, dist_tot, xi, S_vert);
    

    return S_vert;
} 

// void interp_ray_positions(rayF** raylist, double n_x, double n_y, double n_z, int t_ind, rayT* rayout) {
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

    // print_vector(vector<double>(W, W+8));
    // print_array(W, 8);
    for (int jj=0; jj<8; jj++){  // Corner rays
    // print_array(raylist[jj]->pos[t_ind].data(), 3);
        // Vector-valued
        for (int ii=0; ii<3; ii++){  // X, Y, Z
            // (rayout->pos)[ii]   += W[jj]*((raylist[jj]->pos[t_ind]).data()[ii]);
            // (rayout->vgrel)[ii] += W[jj]*((raylist[jj]->pos[t_ind]).data()[ii]);
            (rayout->n)[ii]     += W[jj]*((framelist[jj].n).data()[ii]);
            (rayout->B0)[ii]    += W[jj]*((framelist[jj].B0).data()[ii]);

        }

        // scalar-valued here
        // cout << "corner w: " << raylist[jj]->w << "\n";
        rayout->time +=  W[jj]*(framelist[jj].time);
        rayout->w += W[jj]*(framelist[jj].w);
        rayout->inp_pwr += W[jj]*(framelist[jj].inp_pwr);
        rayout->damping += W[jj]*(framelist[jj].damping);

        rayout->stixR += W[jj]*(framelist[jj].stixR);
        rayout->stixL += W[jj]*(framelist[jj].stixL);
        rayout->stixP += W[jj]*(framelist[jj].stixP);
        rayout->stixS += W[jj]*(framelist[jj].stixS);
        rayout->stixD += W[jj]*(framelist[jj].stixD);
        // rayout->stixA += W[jj]*(framelist[jj].stixA);
        // rayout->stixB += W[jj]*(framelist[jj].stixB);


    }


    // cout << "in: ";
    // print_array(rayout->pos, 3);
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


        // k = n_vec*w/C;
        // kmag = k.norm();
        // Bhat = B0.array()/B0.norm();
        // kpar = k.dot(Bhat); //k.array()*Bhat.array();
        // kperp = (k - kpar*Bhat).norm();

        // // Theta is the angle between parallel and perpendicular K
        // theta = atan2(-kperp, -kpar);   // negation on both sides matches the original raytracer
        //                                 // (Probably means the magnitude of the b-field is backwards
        //                                 // between forrest or jacob)
        // // Some trig.
        // sin_th = sin(theta);
        // cos_th = cos(theta);
        // sin_th_sq = pow(sin_th,2);
        // cos_th_sq = pow(cos_th,2);

        // A = S*sin_th_sq + P*cos_th_sq;
        // B = R*L*sin_th_sq + P*S*(1+cos_th_sq);



        ray->stixR.push_back(R);
        ray->stixL.push_back(L);
        ray->stixP.push_back(P);
        ray->stixS.push_back(S);
        ray->stixD.push_back(D); 
        // ray->stixA.push_back(A);
        // ray->stixB.push_back(B);

        // --------------------------------------------

    } // ii (timestep)
}


void init_EA_array(EA_segment* EA_array, double lat, double lon, int itime_in[2], int model_number) {

    double x_in[3], x_in_geocar[3], x_sm[3];         
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



    x_in = {1, lat, lon};

    // Get start coordinate in SM cartesian:
    x_in_geocar = {x_in[0], x_in[1], x_in[2]};
    degcar(x_in_geocar);
    mag_to_sm_d_(itime_in, x_in_geocar, x_sm);

    cout << "orig: ";
    print_array(x_in, 3);
    cout << "SM: ";
    print_array(x_sm, 3);

    
    // This works!
    //Nsteps = trace_fieldline(x_in, x_fl, TRACER_STEP);

    double tsyg_params[10] = {0};
    double VG[3];
   
    // int model_number = 0;

     // Setup for IGRF:    
    load_TS05_params(itime_in, tsyg_params, VG);
    init_igrf(itime_in);

    // set DST:
    tsyg_params[1] = -20;

    // Trace field line:
    Nsteps = trace_fieldline(itime_in, x_sm, x_fl, b_fl, TRACER_STEP, model_number, tsyg_params);
    

    // FILE* dump;

    // dump = fopen("trace.txt", "w");

    // for (int i=0; i < Nsteps; i++ ) {
    //     fprintf(dump, "%i %g %g %g %g\n",i, x_fl[i][0], x_fl[i][1], x_fl[i][2], b_fl[i]);
    // }
    // fclose(dump);


    double lats[Nsteps];
    double dist_n[Nsteps];
    double ftc_n[Nsteps];   // Flight time constant to northern hemi

    double target_lat = EALimN;
    int EA_index = 0;
    
    // Get magnetic latitudes of entries
    for (int i=0; i < Nsteps; i++) {
        sm_to_mag_d_(itime_in, x_fl[i], x_cur);
        cardeg(x_cur);
        lats[i] = x_cur[1];

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

    cout << "total distance: " << dist_n[Nsteps] << "\n";
    // Get effective L-shell
    // (Do we take this to be the radius at geomag equator, or the maximum?)
    int Lsh_index = nearest(lats, Nsteps, 0, true);

    double Lsh = norm(x_fl[Lsh_index], 3);
    cout << "L shell: " << Lsh << "\n";

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
    cout << "Alpha_eq: " << alpha_eq << "\n";

    // Calculate flight-time constants:
    int iono_ind = 0;
    double alt = 0;
    // Find index of ionosphere
    while (alt < (H_IONO/R_E + 1)) {
        iono_ind +=1;
        Vector3d xtmp = Map<VectorXd>(x_fl[iono_ind],3,1);
        alt = xtmp.norm();
    }

    cout << "Iono ind: " << iono_ind << "\n";

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



    double targ_lat = EALimN;

    // generate entries for each EA segment
    for (int i=0; i < NUM_EA; i++) {
        ind = nearest(lats, Nsteps, targ_lat, true);
        cout << "lat: " << lats[ind] << "\n";
        cout << "Lsh: " << Lsh << "\n";
        // L shell:
        EA_array[i].Lsh = Lsh;

        // Latitude:
        EA_array[i].lat = lats[ind];

        // Distance to northern and southern ionosphere intersections:
        EA_array[i].dist_to_n = dist_n[ind]                    - H_IONO/R_E;
        EA_array[i].dist_to_s = dist_n[Nsteps-1] - dist_n[ind] - H_IONO/R_E;

        // Flight-time constants:
        EA_array[i].ftc_n     = ftc_n[ind];
        EA_array[i].ftc_s     = ftc_n[Nsteps-1] - ftc_n[ind];
        cout << "ftc n: " << EA_array[i].ftc_n << " ftc s: " << EA_array[i].ftc_s << "\n";

        // Intersection point:
        EA_array[i].ea_pos = Map<VectorXd>(x_fl[ind],3,1);
        // Get B:
        // Bomag = b_fl[ind];
        // cout << "Bo: " << Bomag << "\n";
        bmodel(itime_in, x_fl[ind], tsyg_params, model_number, Bo);

        // Get unit vector pointing along field line:
        // (i.e., normal vector to this EA segment)
        Bomag = norm(Bo, 3);
        // Bomag = b_fl[ind];
        EA_array[i].ea_norm = Map<VectorXd>(Bo, 3, 1)/Bomag;

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
        
        EA_array[i].radius = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))/2.0; 
        printf("EA_length: %g\n",2*R_E*EA_array[i].radius);


        // Calculate electron gyrofrequency:
        EA_array[i].wh = (Q_EL/M_EL)*Bomag;

        // Spatial derivative of gyrofrequency:
        // (Currently evaluating using Jacob's formula -- dipole model assumption!!)
        slat_term = sqrt(1+3*slam*slam);

        EA_array[i].dwh_ds = 3*(EA_array[i].wh)/(Lsh*R_E)*slam/slat_term*
                (1/(slat_term*slat_term) + 2/(clam*clam));



        // Calculate local loss cone:
        EA_array[i].alpha_lc = asin(sqrt(Bomag/norm(B_iono,3)));
        
        // Save equatorial loss cone:
        EA_array[i].alpha_eq = alpha_eq;

        // dv_|| / ds:  (again, Jacob's formula -- dipole model assumption)
        EA_array[i].dv_para_ds = -0.5*pow(sin(EA_array[i].alpha_lc),2)/
                                      cos(EA_array[i].alpha_lc)/
                                      EA_array[i].wh*EA_array[i].dwh_ds;


        // // analytical (for comparison)
        // double epsm = (1/Lsh)*(R_E+H_IONO)/R_E;
        // double alpha_eq = asin(sqrt( pow(epsm,3)/sqrt(1+3*(1-epsm)) ));
        // double slat_term = sqrt(1+3*slam*slam);
        // double alpha_calc = asin(sqrt( slat_term/pow(clam,6) )*sin(alpha_eq));

        // cout << "alpha_lc: " << EA_array[i].alpha_lc << " calculated: " << alpha_calc << "\n";


        // Bump index
        targ_lat-= EAIncr;
    }


    // Loop through again to calculate ds:
    // (Distance along field line between EA segments)
    for (int i=0; i < NUM_EA; i++) {
        if (i==0) {
            EA_array[i].ds = EA_array[i+1].dist_to_n - EA_array[i].dist_to_n;
        } else if (i== NUM_EA - 1) {
            EA_array[i].ds = EA_array[i].dist_to_n - EA_array[i-1].dist_to_n;
        } else {
            EA_array[i].ds = (EA_array[i+1].dist_to_n - EA_array[i].dist_to_n)*0.5
                            +(EA_array[i].dist_to_n - EA_array[i-1].dist_to_n)*0.5;
        }

        EA_array[i].ds *= R_E;  // Meters

        // clam = cos(EA_array[i].lat*D2R);
        // slam = sin(EA_array[i].lat*D2R);
        // slat_term = sqrt(1+3*slam*slam);
        // double ds_calc = EA_array[i].Lsh*slat_term*clam*EAIncr*D2R*R_E; 

        // cout << "lat: " << EA_array[i].lat << " ds: " << EA_array[i].ds << " calc: " << ds_calc << 
        //      "dist to n: " << EA_array[i].dist_to_n << "\n";
    }
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
            lside = longitude_interval(r1[2], EAr[2]);
        }
        sl1 = longitude_interval(r1[2], EAr[2]);
        sl2 = longitude_interval(r2[2], EAr[2]);

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


double longitude_interval(double ra, double r0) {
    // Pretty sure this won't wrap around happily, but it works for now.
    // (10.28.16)

    // cout << "ra: " << ra*R2D << " r0: " << r0*R2D << " diff: " << R2D*(ra - r0) << "\n";

    const double width = 2.5;

    if ( R2D*(ra - r0) < - width) {
        return -1;
    } else if ( R2D*(ra - r0) > width) {
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


void dump_EA_array(EA_segment EA_array[NUM_EA], string filename) {
    // Write the EA array to a file, so I can plot it.

    FILE * file;

        // Save it
    cout << "saving to file " << filename << "\n";
    file = fopen(filename.c_str(), "w");

    if (file != NULL) {
        
        for (int i=0; i < NUM_EA; i++) {
            fprintf(file, "%g %g %g %g %g %g %g\n",
                EA_array[i].ea_pos[0], EA_array[i].ea_pos[1],  EA_array[i].ea_pos[2],
                EA_array[i].ea_norm[0],EA_array[i].ea_norm[1], EA_array[i].ea_norm[2],
                EA_array[i].radius);
        }
    } else {
        cout << "Could not open file " << filename.c_str() << "\n";
    }
}


void calc_resonance(rayT* ray, EA_segment* EA, double v_tot_arr[NUM_E], 
    double da_N[NUM_E][NUM_TIMES], double da_S[NUM_E][NUM_TIMES]) {
    // Performs resonance calculation for a single ray entry, and records
    // changes in pitch angle into da_N and da_S (for northern and southern hemis)
    double t, w, pwr;
    Vector3d n, B0;
    double psi, mu, mu_sq, spsi, cpsi, spsi_sq, cpsi_sq;
    double k, kx, kz;
    double n_x, n_z;
    double wh, alpha_lc, calph, salph, ds, dv_para_ds, dwh_ds;
    double Y;
    double stixS, stixD, stixA, stixB, stixX, stixR, stixL, stixP;
    double rho1, rho2, Byw_sq;
    double Byw, Exw, Eyw, Ezw, Bxw, Bzw;
    double R1, R2, w1, w2, alpha1;
    double t1, t2, t3;
    double direction;
    double v_para_res, v_tot_res, E_res;
    double e_starti, e_endi;
    double slat, clat, slat_term;
    double alpha_eq;
    double ftc_n, ftc_s;
    
    double v_tot, v_para, v_perp;
    double gamma, alpha2, beta, wtau_sq, T1;
    double eta_dot, dalpha_eq;
    double v_para_star, v_para_star_sq;
    double AA, BB;
    double Farg, Farg0, Fs, Fs0, Fc, Fc0;
    double dFs_sq, dFc_sq;
    double dalpha, alpha_eq_p;
    double flt_time;

    Vector3d kvec, Bhat;
    int timei;

    // Parameters from the EA array (maybe we should calculate them here, for organization)
    wh = EA->wh;
    alpha_lc = EA->alpha_lc;
    calph = cos(alpha_lc);
    salph = sin(alpha_lc);
    ds    = EA->ds;
    dv_para_ds = EA->dv_para_ds;
    dwh_ds     = EA->dwh_ds;
    slat  = sin(D2R*EA->lat);
    clat  = cos(D2R*EA->lat);
    alpha_eq = EA->alpha_eq;
    ftc_n = EA->ftc_n;
    ftc_s = EA->ftc_s;


    t = ray->time + ray->dt/2.;
    w = ray->w + FREQ_STEP_SIZE*PI;
    pwr = (ray->inp_pwr)*(ray->damping)/(1 + ray->num_rays)*FREQ_STEP_SIZE;

    n  = ray->n /(1 + ray->num_rays);
    B0 = ray->B0;


    if (pwr < WAVE_PWR_THRESH) {
        return;
    }

    // For this calculation, we're working in a frame with z parallel to 
    // the background magnetic field. 

    // Angle between wavenormal and background magnetic field.
    // ( -90 - (dot product) matches original raytracer's output... 11.1.2016)

    // psi = -90*D2R - n.dot(B0)/(n.norm()*B0.norm());


    kvec = n*w/C;
    k    = kvec.norm();
    Bhat = B0.array()/B0.norm();
    kz = -1*kvec.dot(Bhat); //k.array()*Bhat.array();
    kx = (kvec + kz*Bhat).norm();

    // Theta is the angle between parallel and perpendicular K
    psi = atan2(-kx, kz);



    // printf("t: %g psi: %2.3f  pwr: %g num_rays: %i\n",t, psi*R2D, pwr, ray->num_rays);
    // cout << "t: " << t << " psi: " << psi*R2D << " Num_rays: "<< ray->num_rays <<  "\n";



    slat_term = sqrt(1+3*slat*slat);

    // cout << "wh: " << wh << "\n";
    // cout << "t= " << t << " psi= " << R2D*psi << "\n";
    mu = n.norm();
    mu_sq = pow(mu, 2);
    spsi = sin(psi);
    cpsi = cos(psi);
    spsi_sq = pow(spsi, 2);
    cpsi_sq = pow(cpsi, 2);

    n_x =  mu*fabs(spsi);
    n_z =  mu*cpsi;

    // k = w*mu/C;
    // kx = w*n_x/C;
    // kz = w*n_z/C;

    Y = wh / w;


    // printf("kx: %g kperp: %g kz: %g kpar: %g dkx: %g dkz: %g\n",kx, kperp, kz, kpar, kx - kperp, kz + kpar);

    stixP = ray->stixP /(1 + ray->num_rays);
    stixR = ray->stixR /(1 + ray->num_rays);
    stixL = ray->stixL /(1 + ray->num_rays);
    stixS = ray->stixS /(1 + ray->num_rays);
    stixD = ray->stixD /(1 + ray->num_rays);



    stixA = stixS*spsi_sq       + stixP*cpsi_sq;
    stixB = stixR*stixL*spsi_sq + stixP*stixS*(1+cpsi_sq);

    stixX = stixP/(stixP - mu_sq*spsi_sq);  // Ristic 3.2, pg 41

    // Polarization ratios
    rho1 = ((mu_sq-stixS)*mu_sq*spsi*cpsi)/(stixD*(mu_sq*spsi_sq-stixP));
    rho2 = (mu_sq - stixS) / stixD ;

    // (bortnik 2.28)
    Byw_sq =  (2.0*MU0/C) * pwr * stixX*stixX * rho2*rho2 * mu*fabs(cpsi)/
                       sqrt(  pow((tan(psi)-rho1*rho2*stixX),2) + 
                              pow( (1+rho2*rho2*stixX), 2 ) );

    // RMS wave components
    Byw = sqrt(Byw_sq);
    Exw = fabs(C*Byw * (stixP - n_x*n_x)/(stixP*n_z)); 
    Eyw = fabs(Exw * stixD/(stixS-mu_sq));
    Ezw = fabs(Exw *n_x*n_z / (n_x*n_x - stixP));
    Bxw = fabs(Exw *stixD*n_z /C/ (stixS - mu_sq));
    Bzw = fabs((Exw *stixD *n_x) /(C*(stixX - mu_sq)));
    
    // printf("Byw: %g Exw: %g Eyw: %g Ezw: %g Bxw: %g Bzw: %g\n",
    //           Byw,    Exw,    Eyw,    Ezw,    Bxw,    Bzw);
    // cout << " Byw_sq: " << Byw_sq << " pwr: " << pwr << " damping: " << ray->damping << " inp: " << ray->inp_pwr << "\n";
    // // printf("\nByw_sq: %g, rho1: %g, rho2: %g, \nstixS: %g, stixB: %g\n", Byw_sq, rho1, rho2, stixS, stixB);  
    // printf("\nn_x: %g, n_z: %g, stixX: %g, stixD: %g, stixA: %g,mu_sq: %g\n", n_x, n_z, stixX, stixD, stixA, mu_sq);
      
    // Oblique integration quantities
    R1 = (Exw + Eyw)/(Bxw+Byw);
    R2 = (Exw - Eyw)/(Bxw-Byw);
    w1 = Q_EL/(2*M_EL)*(Bxw+Byw);
    w2 = Q_EL/(2*M_EL)*(Bxw-Byw);
    alpha1 = w2/w1;

    // Sum over all resonance modes
    for (int mres=-SCATTERING_RES_MODES; mres <= SCATTERING_RES_MODES; mres++) {
        // get parallel resonance velocity
        t1 = w*w*kz*kz;
        t2 = pow((mres*wh),2)-w*w;
        t3 = kz*kz + pow((mres*wh),2)/(pow(C*cos(alpha_lc),2));
        if(mres==0) {
          direction = -1*signof(kz); //kz/fabs(kz);
        } else {
          direction = signof(kz)*signof(mres); //kz/fabs(kz) * mres/fabs(mres) ;
        }

        v_para_res = ( direction*sqrt(t1 + t2*t3) - w*kz ) / t3;
        v_tot_res = v_para_res / cos(alpha_lc); 
        E_res = E_EL*( 1.0/sqrt( 1.0-(v_tot_res*v_tot_res/(C*C)) ) -1.0 );

        // indexes into the energy / velocity arrays
        e_starti = floor((log10(E_res) - E_EXP_BOT - 0.3)/(DE_EXP));
        e_endi   =  ceil((log10(E_res) - E_EXP_BOT + 0.3)/(DE_EXP));

        if(e_endi>NUM_E) e_endi=NUM_E;
        if(e_starti>NUM_E) e_starti=NUM_E;
        if(e_endi<0) e_endi=0;
        if(e_starti<0) e_starti=0;


        // begin V_TOT loop here
        for(int e_toti=e_starti; e_toti < e_endi; e_toti++) {
            v_tot = direction*v_tot_arr[e_toti];
            v_para = v_tot * calph;
            v_perp = fabs(v_tot * salph);

            gamma = 1.0 / sqrt(1 - pow((v_tot/C),2));   // Relativisitc correction
            alpha2 = Q_EL*Ezw /(M_EL*gamma*w1*v_perp);
            beta = kx*v_perp / wh ;
            wtau_sq = pow((-1),(mres-1)) * w1/gamma * 
                ( jn( (mres-1), beta ) - 
                alpha1*jn( (mres+1) , beta ) +
                gamma*alpha2*jn( mres , beta ) ); 
            T1 = -wtau_sq*(1+ ( (calph*calph) / (mres*Y-1) )  );

            // Now -- start the analytical evaluation!!
            if( fabs(EA->lat) < 1e-3) {
                // Near the equator we can use a simplified expression:

                eta_dot = mres*wh/gamma - w - kz*v_para;

                if(fabs(eta_dot)<10) {
                    // Bortnik A.31
                    dalpha_eq = fabs(T1/v_para)*ds/sqrt(2);
                } else {
                    // Bortnik A.30
                    dalpha_eq = fabs(T1/eta_dot)*sqrt(1-cos(ds*eta_dot/v_para)); 
                }

            } else { // Full analytical expression required:

                v_para_star = v_para - dv_para_ds*ds/2.0;   // Bortnik A.17
                v_para_star_sq = v_para_star * v_para_star;

                // Bortnik A.18 -- part A1
                AA = - mres/(2.0*v_para_star_sq*gamma)*wh*dv_para_ds
                     + (mres/(2.0*v_para_star*gamma))*dwh_ds * (1 + ds/(2.0*v_para_star)*dv_para_ds) 
                     + w/(2.0*v_para_star_sq)*dv_para_ds;

                // // Bortnik A.18 -- part A0   -- THIS DOES NOT MATCH THE THESIS
                // BB =   mres/(gamma*v_para_star)*wh 
                //      - mres/(gamma*v_para_star)*dwh_ds*(ds/2.0)
                //      - w/v_para_star - kz;

                // Bortnik A.18 -- part A0
                BB =   mres*wh/(gamma*v_para_star)
                     - mres/(gamma*v_para_star)*dwh_ds*(ds/2.0) * (w/v_para_star)*kz;


                // Evaluate Bortnik A.26 -- integration performed thru Fresnel functions
                Farg = (BB + 2*AA*ds) / sqrt(2*PI*fabs(AA));
                Farg0 = BB / sqrt(2*PI*fabs(AA));  

                Fresnel(Farg,  &Fs,  &Fc);
                Fresnel(Farg0, &Fs0, &Fc0);

                dFs_sq = pow((Fs - Fs0),2);
                dFc_sq = pow((Fc - Fc0),2);

                dalpha = sqrt(PI/4/fabs(AA))*fabs(T1/v_para)*sqrt(dFs_sq+dFc_sq);

                // Map the local change in pitch angle to the equivalent
                // pitch angle at the equator:  (still using dipole model here...)
                alpha_eq_p = asin( sin(alpha_lc+dalpha)*pow(clat,3) / 
                                   sqrt(slat_term) );
                dalpha_eq = alpha_eq_p - alpha_eq;

                // Determine where to bin the output:
                if(direction>0) {
                    flt_time = fabs(ftc_n/v_para);
                } else {
                    flt_time = fabs(ftc_s/v_para);
                }
                

                // Get time index into output array
                timei = round((t + flt_time)/TIME_STEP);
                
                // cout << "flt_time: " << flt_time << " t: " << t << " timei: " << timei << "\n";
                // Save it!
                if (direction > 0) {
                    da_N[e_toti][timei] += dalpha_eq*dalpha_eq;
                } else {
                    da_S[e_toti][timei] += dalpha_eq*dalpha_eq;
                }
            }   
        }   // v_tot, e_tot
    } // Resonant modes

}


void add_rayT(rayT* rayA, rayT* rayB) {

    rayA->inp_pwr += rayB->inp_pwr;       
    rayA->damping += rayB->damping;      
    rayA->n       += rayB->n;         
    rayA->B0      += rayB->B0;        
    rayA->stixP   += rayB->stixP;         
    rayA->stixR   += rayB->stixR;         
    rayA->stixL   += rayB->stixL;         
    // rayA->stixA   += rayB->stixA;         
    // rayA->stixB   += rayB->stixB;         
    rayA->stixS   += rayB->stixS;         
    rayA->stixD   += rayB->stixD;

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
    // cout << "Scalars\n";
    // Scalar-valued
    frame->damping = ( rayfile->damping[iMid+1]-rayfile->damping[iMid] )*M + rayfile->damping[iMid];
    frame->stixP = ( rayfile->stixP[iMid+1]-rayfile->stixP[iMid] )*M + rayfile->stixP[iMid];
    frame->stixR = ( rayfile->stixR[iMid+1]-rayfile->stixR[iMid] )*M + rayfile->stixR[iMid];
    frame->stixL = ( rayfile->stixL[iMid+1]-rayfile->stixL[iMid] )*M + rayfile->stixL[iMid];
    frame->stixS = ( rayfile->stixS[iMid+1]-rayfile->stixS[iMid] )*M + rayfile->stixS[iMid]; 
    frame->stixD = ( rayfile->stixD[iMid+1]-rayfile->stixD[iMid] )*M + rayfile->stixD[iMid];
    // frame->stixA = ( rayfile->stixA[iMid+1]-rayfile->stixA[iMid] )*M + rayfile->stixA[iMid];
    // frame->stixB = ( rayfile->stixB[iMid+1]-rayfile->stixB[iMid] )*M + rayfile->stixB[iMid];   
    

    // Stuff that doesn't need interpolation:
    frame->w    = rayfile->w;
    frame->time = t_target;
    frame->inp_pwr = rayfile->inp_pwr;
    // cout << t_target << " " << frame->damping <<"\n";

}


