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
    // cout << " gc_dist: " << gc_distance << "\n";
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


void interp_ray_fine(rayF** raylist, double n_x, double n_y, double n_z, int t_ind, rayT* out) {
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

    for (int jj=0; jj<8; jj++){  // Corner rays

        // Vector-valued
        for (int ii=0; ii<3; ii++){  // X, Y, Z
            out->pos[ii] += W[jj] * (raylist[jj]->pos[t_ind][ii]);
            out->vgrel[ii] += W[jj] * (raylist[jj]->pos[t_ind][ii]);
        }
        // scalar-valued here
        // cout << "corner w: " << raylist[jj]->w << "\n";
        out->w += W[jj]*(raylist[jj]->w);

        out->stixR += W[jj]*(raylist[jj]->stixR[t_ind]);
        out->stixL += W[jj]*(raylist[jj]->stixL[t_ind]);
        out->stixP += W[jj]*(raylist[jj]->stixP[t_ind]);
        out->stixS += W[jj]*(raylist[jj]->stixS[t_ind]);
        out->stixD += W[jj]*(raylist[jj]->stixD[t_ind]);
        out->stixA += W[jj]*(raylist[jj]->stixA[t_ind]);
        out->stixB += W[jj]*(raylist[jj]->stixB[t_ind]);
    }

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
    double EA_a, EA_b, EA_c;


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
    }

    cout << "total distance: " << dist_n[Nsteps] << "\n";
    // Get effective L-shell
    //  (Do we take this to be the radius at geomag equator, or the maximum?)
    int Lsh_index = nearest(lats, Nsteps, 0, true);

    double Lsh = norm(x_fl[Lsh_index], 3);
    cout << "L shell: " << Lsh << "\n";

    int ind;
    double targ_lat = EALimN;


    // generate entries for each EA segment
    for (int i=0; i < NUM_EA; i++) {
        ind = nearest(lats, Nsteps, targ_lat, true);
        cout << "lat: " << lats[ind] << "\n";

        // L shell:
        EA_array[i].Lsh = Lsh;

        // Distance to northern and southern ionosphere intersections:
        EA_array[i].dist_to_n = dist_n[i]                  - (R_E + H_IONO)/R_E;
        EA_array[i].dist_to_n = dist_n[Nsteps] - dist_n[i] - (R_E + H_IONO)/R_E;

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

        // Bump index
        targ_lat-= EAIncr;
    }


    // double x_tmp[3];
    // // Read out found data:
    // for (int j=0; j < NUM_EA; j++) {
    //     sm_to_mag_d_(itime_in, EA_array[j].ea_pos.data(), x_tmp);
    //     cardeg(x_tmp);
    //     print_array(x_tmp, 3);

    //     cout << EA_array[j].radius << "\n";

    // }

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
    Vector3d p0;
    Vector3d n;
    Vector3d l0, l1;

    double r0, r1, r2;

    p0 = EA.ea_pos;
    n  = EA.ea_norm;
    r0 = EA.radius;

    // // Check each of the 8 rays:
    for (int rr = 0; rr < 8; rr++) {

        l0 = Map<VectorXd>(cur_rays[rr]->pos[t - 1].data(), 3, 1);
        l1 = Map<VectorXd>(cur_rays[rr]->pos[t].data(), 3, 1);

        r1 = (l0 - p0).norm();
        r2 = (l1 - p0).norm();

        // cout << "r1: " << r1 << " r2: " << r2 << " EA: " << EA.radius << "\n";
        if ( (r1 < r0) ) {
            cout << "r1: " << r1 << " r2: " << r2 << " EA: " << EA.radius << "\n";
            return true;
        }
    }

    return false;
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


