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

    double x_in[3];         
    double x_out[3];

    double b_out[3];

    int itime_in[2];
    itime_in[0] = 1000*iyr + idoy;
    itime_in[1] = isec*1e3;


    x_in = {1, lat, lon};

    cout << "orig: ";
    print_array(x_in, 3);

    // double r = 1;
    // double latin = lat*D2R;
    // double lonin = lon*D2R;

    // pol_to_cart_d_(&latin, &lonin, &r, x_out);

    // cout << "xformd: ";
    // print_array(x_out, 3);

    // sph_car_(&r, &lat, &lon, x_out);


    // cout << "IRBEM: ";
    // print_array(x_out, 3);

    // degcar(x_in);
    // cout << "cart: ";
    // print_array(x_in, 3);

    // cardeg(x_in);
    // cout << "back: ";
    // print_array(x_in, 3);

    int Nsteps;
    // bmodel_dipole(x_in, b_out);
    // cout << "b (sph): ";
    // print_array(b_out, 3);

    // transform_data_sphcar(b_out, lat, lon);



    double x_fl[TRACER_MAX][3];
    double ds = 0.0005;
    Nsteps = trace_fieldline(x_in, x_fl, ds);

    for (int i=0; i < Nsteps; i++) {
        cout << "x(" << i << ") : ";
        print_array(x_fl[i],3);
    }

    // cout << "B: ";
    // print_array(b_out, 3);

    // // Loop through possible EA segments:
    // for (int ii=0; ii < NUM_EA; ii++) {

    // }


}

