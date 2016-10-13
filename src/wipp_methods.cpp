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


void init_EA_array(EA_segment* EA_array, double lat, double lon, int iyr, int idoy, double isec) {

    // First: trace the field line from output location:
    double kext = 0;    // External field model (0 for none)
    int options[5];
        options[0] = 1;     // Compute L* and phi
        options[1] = 30;     // Interval in days between updating the IGRF model (0 = 1 year)
        options[2] = 0;     // Resolution to compute L* to [0..9] (0 suggested)
        options[3] = 0;     // Another resolution setting (0 suggested)
        options[4] = 0;     // Internal mag. field model (0 = igrf, 1 = eccentric tilted dipole)
    int sysaxes = 7;       // Coordinate system of input data (6 = MAG cartesian)
    double ds = 1e-3;       // Integration step size along field line (in Earth radii)

    double maginput[25] = {0};    // Model parameters

    double posit[3000][3];  // Output space -- up to 3000 vector coordinates
    int Nposit;          // Actual number of entries

    double x_in[3];         // inputs to the field line tracer (mag cartesian, Re)
    double x_out[3];

    int itime_in[2];
    itime_in[0] = 1000*iyr + idoy;
    itime_in[1] = isec*1e3;

    double radius = 1; //(R_E + H_IONO)/R_E;

    double R0 = 1;

    double olat, olon, orad;
    double lm;
    double blocal[3000];
    double bmin;
    double xj;
    double lat_in, lon_in;

    // lat_in = D2R*lat;
    // lon_in = D2R*lon;
    // pol_to_cart_d_(&lat_in, &lon_in, &radius, x_in);

    // trace_field_line_towards_earth1_(&kext, options, &sysaxes, 
    //             &iyr, &idoy, &isec,
    //             x_in, x_in + 1, x_in + 2,
    //             maginput, &ds, posit, &Nposit);

    x_in[0] = 1;
    x_in[1] = 45;
    x_in[2] = 0;
    

    // sph_car_(&radius, &lat, &lon, x_in);

    cout << "calling field line duder\n";
    trace_field_line2_1_(&kext, options, &sysaxes,
                &iyr, &idoy, &isec,
                x_in, x_in + 1, x_in + 2,
                maginput, &R0, &lm, 
                blocal, &bmin, &xj,
                posit, &Nposit);



    cout << "recieved " << Nposit << " elements\n";
    // Write the output for debugging:

    cout << "R0: " << R0 << " Lm: " << lm;
    cout << " blocal: " << blocal[0] << " bmin: " << bmin << " xj: " << xj << "\n";


        cout << "Start: " << posit[0][0] << ", " << posit[0][1] << ", " << posit[0][2] << "\n";
        cout << "End: "   << posit[Nposit-1][0] << ", " << posit[Nposit - 1][1] << ", " << posit[Nposit -1][2] << "\n"; 

             // IRBEM cast
            geo2mag1_(&iyr, posit[0], x_out);
            car_sph_(x_out, &orad, &olat, &olon);

            // olat = R2D*olat; olon = R2D*olon;
            cout << "Start: " << orad << ", " << olat << ", " << olon <<"\n";

             // IRBEM cast
            geo2mag1_(&iyr, posit[Nposit - 1], x_out);
            car_sph_(x_out, &orad, &olat, &olon);
            cout << "End: " << orad << ", " << olat << ", " << olon <<"\n";

            // // olat = R2D*olat; olon = R2D*olon;
            // cout << "End: " << orad << ", " << olat << ", " << olon <<"\n";



    // FILE *file = fopen("/shared/users/asousa/WIPP/3dWIPP/outputs/fieldlinelog.txt", "w");    

    // if (file != NULL) {
    //     cout << "logging\n";

        for (long i=0; i < Nposit; i++) {

        cout << posit[i][0] << ", " << posit[i][1] << ", " << posit[i][2] << "\n";

            // Xformd cast
            // geo_to_mag_d_(itime_in, posit[i], x_out);
            // cart_to_pol_d_(x_out, &olat, &olon, &orad);

            // // IRBEM cast
            // geo2mag1_(&iyr, posit[i], x_out);
            // car_sph_(x_out, &orad, &olat, &olon);

            // // olat = R2D*olat; olon = R2D*olon;
            // cout << "i: " << i << " -- " << orad << ", " << olat << ", " << olon <<"\n";
            // // cout << "i: " << i << " -- " << posit[i][0] << ", " << posit[i][1] << ", " << posit[i][2] <<"\n";

        }

    // } else {
    //     cout << "something's fucky\n";
    // }


    // // Loop through possible EA segments:
    // for (int ii=0; ii < NUM_EA; ii++) {

    // }


}

