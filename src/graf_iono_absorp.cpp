#include <wipp.h>

double graf_iono_absorp(float lat, long f, double mlt)
{
    // A model of VLF wave power attenuation between 100 and 1000 km altitude.
    // Based on Graf and Cohen 2013: "Analysis of experimentally validated
    // trans-ionospheric attenuation estimates of VLF signals", figure 7.

    // Data picked from plots and fitted to an exponential function.
    // Interpolates / extrapolates in log space for frequencies other than 2kHz
    // and 20kHz. Day / night variation is weighted using a sigmoid function
    // with a slope arbitrarily set to "looks good".
    // --APS 1.2017 asousa@stanford.edu
    
    const double mltslope = 0.5;

    double mlt_mod = fmod(mlt,24);

    float db2i, db20i;
    double  db2iLog, db20iLog, m, c; 

    double p20D[3] = {215.74122661,   11.9624129,  9.01400095};
    double p20N[3] = {117.54370955,   7.40762459,  0.90050155};
    double p2D[3]  = {55.94274086,    11.91761368, 4.09353494};
    double p2N[3]  = {11.99682851,    9.53682009,  0.23617706};

    double a20D, a20N, a2D, a2N;
    double mD, mN, cD, cN;
    double aD, aN;
    double s1, s2, s;

    
    // # Attenuation values at 20kHz and 2kHz, day and night
    a20D  = log10(p20D[0]* exp(-fabs(lat)/p20D[1]) + p20D[2]);
    a20N  = log10(p20N[0]* exp(-fabs(lat)/p20N[1]) + p20N[2]);
    a2D   = log10(p2D[0] * exp(-fabs(lat)/p2D[1]) + p2D[2]);
    a2N   = log10(p2N[0] * exp(-fabs(lat)/p2N[1]) + p2N[2]);   
  
    // Interpolate / extrapolate between 20kHz and 2kHz values (log space)
    // Should be m'=m/(log10(20)-log10(2)) but denominator = 1
    mD = (a20D - a2D); // day slope
    mN = (a20N - a2N); // night slope
    cD = (a2D + a20D)/2. - mD*0.8010; // day offset
    cN = (a2N + a20N)/2. - mN*0.8010; // night offset
    // (this is just (log10(2)+log10(20))/2 )

    // now extrapolate to desired f value in log10 domain
    // get 10^results
    aD = pow(10.0, mD*log10(f/1000.0) + cD);
    aN = pow(10.0, mD*log10(f/1000.0) + cN);

    // # Weight day + night according to MLT (logistic function)
    s1 = 1.0/(1 + exp(((mlt_mod) - 18)/mltslope));
    s2 = 1.0/(1 + exp(((mlt_mod) - 6)/mltslope));
    s = s1 - s2;
    
    // # Select day curve for s = 1, night curve for s = 0
    return s*aD + (1.0 - s)*aN;
}



double total_input_power(double flash_pos_sm[3], double i0, 
                        double latmin, double latmax, double lonmin, double lonmax, double wmin, double wmax, int itime_in[2]) {
    // Determine the total input power tracked by the set of guide rays:

    double tot_pwr = 0;
    double pwr = 0;
    // Integration step sizes
    double dlat = 0.05; 
    double dlon = 0.05;
    double dw   = 5*2*PI;
    double tmp_coords[3] = {0,0,0};
    double x_sm[3];
    double mlt;

    // cout << "now we are here" << endl;
    for (double w = wmin + dw/2; w < wmax; w+=dw) {
        for (double lat = latmin + dlat/2; lat < latmax; lat+=dlat) {
            for (double lon=lonmin; lon < lonmax; lon+=dlon) {
                // cout << "(" << lat << ", " << lon << ")\n";
                tmp_coords = {1 + H_IONO/R_E, lat, lon};
                degcar(tmp_coords);
                mag_to_sm_d_(itime_in, tmp_coords, x_sm);
                mlt = MLT(itime_in, lon);
                pwr = input_power_scaling(flash_pos_sm, x_sm, lat, w, i0, mlt);

                double dist_lat = (R_E + H_IONO)*dlat*D2R;
                double dist_lon = (R_E + H_IONO)*dlon*cos(D2R*lat)*D2R;
                                // Latitude distance      longitude distance       freq dist
                // cout << "dist_lat: " << dist_lat << ", dist_lon: " << dist_lon << "\n";
                tot_pwr += pwr * dist_lat * dist_lon * dw;
            }
        }
    }

    return tot_pwr;
}


double input_power_scaling(double* flash_loc, double* ray_loc, double mag_lat, double w, double i0, double MLT) {
    // Returns ray power at the top of the ionosphere
    // per unit in area and frequency.

    double theta;        // angle between two vectors
    double gc_distance; // great circle distance between two points
    double dist_tot;
    double xi;

    double S, S_vert;
    double attn_factor;
    double w_sq, f;

    double v1[3];
    double v2[3];

    f = w/(2.0*PI);
    
    cardeg(flash_loc, v1);
    cardeg(ray_loc,   v2);
    gc_distance = haversine_distance(v1[1], v1[2], v2[1], v2[2]);

    // total distance up to ionosphere:
    dist_tot = hypot(gc_distance, H_IONO);
    xi = atan2(gc_distance, H_IONO);  // Incident angle

    w_sq =  pow( w , 2 );
    S = ( (1/Z0) * pow( (H_E*i0*2E-7*(sin(xi)/dist_tot)*w*(P_A-P_B)) , 2 ) 
                   /  (  (w_sq+pow(P_A,2))*(w_sq+pow(P_B,2))  )      ) ;
    S_vert = S * cos(xi) ;  // factor for vert prop.

    // Ionosphere absorption model
    attn_factor = pow(10,-(graf_iono_absorp(mag_lat,f, MLT )/10)  );
    S_vert = S_vert * attn_factor;
    // cout << "HEY!" << endl;
    // printf("i0: %2.3f, dist_tot: %2.3f, xi: %2.3f, S_vert; %e\n",i0, dist_tot, xi, S_vert);
    
    return S_vert;
}




