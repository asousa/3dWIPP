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

    double x_fl[TRACER_MAX][3]; // elements along the field line

    double b_dipole[3];    
    int Nsteps;

    int itime_in[2];
    itime_in[0] = 1000*iyr + idoy;
    itime_in[1] = isec*1e3;


    x_in = {1, lat, lon};

    cout << "orig: ";
    print_array(x_in, 3);

    
    // This works!
    //Nsteps = trace_fieldline(x_in, x_fl, TRACER_STEP);


    for (int iyear = iyr; iyear < 2020; iyear++) {

        // int ihr = 0; int imn = 0; int isc = 0;
        // int idoy = 364/2;

        // int ihr = (int)(isec/3600);
        // int imn = int((isec - ihr*3600)/60);
        // int isc = int(isec - ihr*3600 - imn*60);
        // double vgsex = -400.0;
        // double vgsey = 0; double vgsez = 0;


        cout << "iyear: " << iyear << "\n";
        itime_in[0] = 1000*iyear + idoy;

        init_igrf(itime_in);
        igrf_geo(x_in, b_out);


        dipole_geo(itime_in, x_in, b_dipole);



        // recalc_08_(&iyear,&idoy,&ihr,&imn,&isc,
        //          &vgsex, &vgsey, &vgsez);

        // double r = 1;
        // double theta = D2R*(90.0 - lat);
        // double phi   = D2R*(lon);

        // igrf_geo_08_(&r, &theta, &phi,
        //              b_out, b_out + 1, b_out + 2);

        // cout << "b_out: ";
        // print_array(b_out, 3);


        double Bomag = norm(b_out, 3);
        cout << " B_igrf: mag: " << Bomag << " | ";
        print_array(b_out, 3);

        double Bomag_dip = norm(b_dipole, 3);
        cout << " B_dip:  mag: " << Bomag_dip << " | ";
        print_array(b_dipole, 3);


    }


    // Test forward / reverse of data transforms:

    cout << "---------\n";
    
    cout << "orig: "; print_array(b_out,3);

    double x_geo[3];
    
    double btmp[3];
    transform_data_sph2car(x_in[1], x_in[2], b_out, btmp);
    cout << "xyz: "; print_array(btmp,3);

    // transform_data_geo2mag(x_in[1], x_in[2], btmp, b_out);
    // cout << "mag: "; print_array(b_out,3);        

    // transform_data_mag2geo(x_in[1], x_in[2], )

    transform_data_car2sph(x_in[1], x_in[2], btmp, b_out);
    cout << "back: "; print_array(b_out,3);        






//     // Test IGRF model:
//     int ntime = 1;          // Number of points in array
//     int kext = 0;           // External field model
//     int options[5];
//         options[0] = 1;     // Compute L* and Phi
//         options[1] = 30;    // Interval to refresh IGRF, in days
//         options[2] = 0;     // Some precision setting
//         options[3] = 0;     // Another precision setting
//         options[4] = 0;     // Internal field model
//     int sysaxes = 0;        // MAG

//     double maginput[25] = {0};

//     // degcar(x_in);
//     double Bgeo[3];
//     double Bl;

//     get_field_multi_(&ntime, &kext, options, &sysaxes,
//                     &iyr, &idoy, &isec,
//                     x_in, x_in + 1, x_in + 2,
//                     maginput, Bgeo, &Bl);


//     cout << "Bl: " << Bl <<  "\nBgeo: ";
//     print_array(Bgeo, 3);
//     cout << "\n";

// // -----------------
//     // Let's try some rotations between coordinate frames:
//     double M[3][3];

//     double A1[3] = {1, 0, 0};
//     double A2[3] = {0, 1, 0};
//     double A3[3] = {0, 0, 1};

//     double B1[3];
//     double B2[3];
//     double B3[3];

//     geo_to_mag_d_(itime_in, A1, B1);
//     geo_to_mag_d_(itime_in, A2, B2);
//     geo_to_mag_d_(itime_in, A3, B3);

//     // print_array(B1,3);
//     // print_array(B2,3);
//     // print_array(B3,3);

//     M[0][0] = B1[0]; M[0][1] = B2[0]; M[0][2] = B3[0];
//     M[1][0] = B1[1]; M[1][1] = B2[1]; M[1][2] = B3[1];
//     M[2][0] = B1[2]; M[2][1] = B2[2]; M[2][2] = B3[2];

//     double out[3] = {0};
//     for (int row=0; row < 3; row++) {
//         for (int col=0; col < 3; col++) {
//             out[col] += M[col][row]*Bgeo[row];
//         }
//     }

//     cout << "B (mdip): ";
//     print_array(out, 3);


// // -------------------





// // -----------------
//     // Let's try some rotations between coordinate frames:
//     double M[3][3];

//     double A1[3] = {1, 0, 0};
//     double A2[3] = {0, 1, 0};
//     double A3[3] = {0, 0, 1};

//     double B1[3];
//     double B2[3];
//     double B3[3];

//     geo_to_mag_d_(itime_in, A1, B1);
//     geo_to_mag_d_(itime_in, A2, B2);
//     geo_to_mag_d_(itime_in, A3, B3);

//     // print_array(B1,3);
//     // print_array(B2,3);
//     // print_array(B3,3);

//     M[0][0] = B1[0]; M[0][1] = B2[0]; M[0][2] = B3[0];
//     M[1][0] = B1[1]; M[1][1] = B2[1]; M[1][2] = B3[1];
//     M[2][0] = B1[2]; M[2][1] = B2[2]; M[2][2] = B3[2];

//     double out[3] = {0};
//     for (int row=0; row < 3; row++) {
//         for (int col=0; col < 3; col++) {
//             out[col] += M[col][row]*A2[row];
//         }
//     }

//     cout << "B2: ";
//     print_array(B2, 3);
//     cout << "out: ";
//     print_array(out, 3);

// // -------------------



    // for (int i=0; i < Nsteps; i++) {
    //     cout << "x(" << i << ") : ";
    //     print_array(x_fl[i],3);
    // }



}

