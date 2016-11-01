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
    S = ( (1/Z0) * pow( (H_E*i0*2E-7*(sin(xi)/dist_tot)*w*(A-B)) , 2 ) 
                   /  (  (w_sq+pow(A,2))*(w_sq+pow(B,2))  )      ) ;
    S_vert = S * cos(xi) ;  // factor for vert prop.



    attn_factor = pow(10,-(ionoAbsorp(mag_lat,f)/10)  );
    S_vert = S_vert * attn_factor ;

    return S_vert;
} 


void interp_ray_fine(rayF** raylist, double n_x, double n_y, double n_z, int t_ind, rayT* rayout) {
    // raylist: an array of pointers to the 8 corner rays
    // nx, ny, nz: location within the grid (0..1) to interpolate at
    // t_ind: time index



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
            (rayout->pos)[ii]   += W[jj]*((raylist[jj]->pos[t_ind]).data()[ii]);
            (rayout->vgrel)[ii] += W[jj]*((raylist[jj]->pos[t_ind]).data()[ii]);
        }

        // scalar-valued here
        // cout << "corner w: " << raylist[jj]->w << "\n";
        rayout->time +=  W[jj]*(raylist[jj]->time[t_ind]);
        rayout->w += W[jj]*(raylist[jj]->w);
        rayout->inp_pwr += W[jj]*(raylist[jj]->inp_pwr);
        rayout->damping += W[jj]*(raylist[jj]->damping[t_ind]);

        rayout->stixR += W[jj]*(raylist[jj]->stixR[t_ind]);
        rayout->stixL += W[jj]*(raylist[jj]->stixL[t_ind]);
        rayout->stixP += W[jj]*(raylist[jj]->stixP[t_ind]);
        rayout->stixS += W[jj]*(raylist[jj]->stixS[t_ind]);
        rayout->stixD += W[jj]*(raylist[jj]->stixD[t_ind]);
        rayout->stixA += W[jj]*(raylist[jj]->stixA[t_ind]);
        rayout->stixB += W[jj]*(raylist[jj]->stixB[t_ind]);


    }


    // cout << "in: ";
    // print_array(rayout->pos, 3);
}


void calc_stix_parameters(rayF* ray) {
    // Calculate Stix parameters for an entire ray (all timesteps)
    // Results are stored in the ray structure.

    double w, whs, wps2;
    double R, L, P, S, D, a, b;
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

        k = n_vec*w/C;
        kmag = k.norm();
        Bhat = B0.array()/B0.norm();
        kpar = k.dot(Bhat); //k.array()*Bhat.array();
        kperp = (k - kpar*Bhat).norm();

        // Theta is the angle between parallel and perpendicular K
        theta = atan2(kperp, kpar);

        // Some trig.
        sin_th = sin(theta);
        cos_th = cos(theta);
        sin_th_sq = pow(sin_th,2);
        cos_th_sq = pow(cos_th,2);


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

        a = S*sin_th_sq + P*cos_th_sq;
        b = R*L*sin_th_sq + P*S*(1+cos_th_sq);




        // if (ii==1) {
        //     cout << "wps2: " << wps2 << " whs: " << whs << "\n";
        //     cout << "Stix Params [0]: ";
        //     cout << R << " ";
        //     cout << L << " ";
        //     cout << P << " ";
        //     cout << S << " ";
        //     cout << D << " ";
        //     cout << a << " ";
        //     cout << b << "\n";
        // }
        // cout << cur_rays[dd]->time.size() << " ";


        ray->stixR.push_back(R);
        ray->stixL.push_back(L);
        ray->stixP.push_back(P);
        ray->stixS.push_back(S);
        ray->stixD.push_back(D);
        ray->stixA.push_back(a);
        ray->stixB.push_back(b);

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

    double b_dipole[3], b_sm[3];    
    int Nsteps;

    double lam, slam, clam, dL_lam, ptR, ptX, ptY, x1, x2, y1, y2;
    double slam2, clam2, rootTerm, x_unit_vect, y_unit_vect, ptL;
    double slat_term;
    double EA_a, EA_b, EA_c;
    double x_iono[3], x_iono_tmp[3];
    double B_iono[3], B_eq[3];
    double alpha_eq;


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
    Nsteps = trace_fieldline(itime_in, x_sm, x_fl, 0.001, model_number, tsyg_params);
    
    double lats[Nsteps];
    double dist_n[Nsteps];

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
    //  (Do we take this to be the radius at geomag equator, or the maximum?)
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

    int ind;
    double targ_lat = EALimN;

    // generate entries for each EA segment
    for (int i=0; i < NUM_EA; i++) {
        ind = nearest(lats, Nsteps, targ_lat, true);
        cout << "lat: " << lats[ind] << "\n";

        // L shell:
        EA_array[i].Lsh = Lsh;

        // Latitude:
        EA_array[i].lat = lats[ind];

        // Distance to northern and southern ionosphere intersections:
        EA_array[i].dist_to_n = dist_n[ind]                  - H_IONO/R_E;
        EA_array[i].dist_to_s = dist_n[Nsteps] - dist_n[ind] - H_IONO/R_E;

        // Intersection point:
        EA_array[i].ea_pos = Map<VectorXd>(x_fl[ind],3,1);

        // Get B:
        bmodel(itime_in, x_fl[ind], tsyg_params, model_number, Bo);

        // Get unit vector pointing along field line:
        // (i.e., normal vector to this EA segment)
        Bomag = norm(Bo, 3);
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
        
        EA_array[i].radius = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))/2; 


        // Calculate electron gyrofrequency:
        EA_array[i].wh = (Q_EL/M_EL)*Bomag;

        // Spatial derivative of gyrofrequency:
        // (Currently evaluating using Jacob's formula -- dipole model assumption!!)
        slat_term = sqrt(1+3*slam*slam);

        EA_array[i].dwh_ds = 3*(EA_array[i].wh)/(Lsh*R_E)*slam/slat_term*
                (1/(slat_term*slat_term) + 2/(clam*clam));



        // Calculate local loss cone:
        EA_array[i].alpha_lc = asin(sqrt(Bomag/norm(B_iono,3)));
        
        // Equatorial loss cone:
        EA_array[i].alpha_eq = alpha_eq;

        // dv_|| / ds: 
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

            lens[i] = trace_fieldline(itime_in, x_sm, grid[i], stepsize, model_number, tsyg_params);
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



bool coarse_mask(rayF* cur_rays[8], int t, EA_segment EA) {
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
        l1 = Map<VectorXd>(cur_rays[rr]->pos[t-1].data(), 3,  1);
        l2 = Map<VectorXd>(cur_rays[rr]->pos[t  ].data(), 3,  1);

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


void calc_resonance(rayT* ray, double v_tot_arr[NUM_E], 
    double da_N[NUM_E][NUM_TIMES], double da_S[NUM_E][NUM_TIMES]) {
    // Performs resonance calculation for a single ray entry, and records
    // changes in pitch angle into da_N and da_S (for northern and southern hemis)

    double t = ray->time + ray->dt/2.;
    double w = ray->w + ray->dw/2.;
    double pwr = (ray->inp_pwr)*(ray->damping);

}



