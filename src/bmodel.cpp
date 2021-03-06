#include <wipp.h>


using namespace std;
using namespace Eigen;


// ---- TILTED DIPOLE MODEL --------------------------------------


void bmodel_dipole(double* x_in, double* B_out) {
    // x_in: R, lat, lon in magnetic dipole coordinates (spherical)
    // B_out: Br, Blat, Blon in magnetic dipole coordinates (Upward, Southward, Eastward)
    // const double Bo = (3.12e-5)*1e9;
    const double Bo = (3.12e-5)*1e9;
    double R, theta, phi;
    double Brad, Btheta, Bphi, Bor3;
    // Walt, pg 30, equations 3.13 - 3.14


    R = x_in[0];
    theta = (90. - x_in[1])*D2R; // Walt defines theta as angle from north (co-latitude)
    phi   = x_in[2]*D2R;

    Bor3 = Bo*pow(R, -3.0);

    Brad = -2.0*Bor3*cos(theta);
    Btheta = -1.0*Bor3*sin(theta);
    Bphi = 0.0;    // Dipole model has no variation in longitude (here for completeness)



    B_out[0] = Brad;        // Up
    B_out[1] = Btheta;      // South
    B_out[2] = Bphi;        // East
}


void dipole_sm(int itime_in[2], double x_in[3], double b_out[3]) {
    //  Dipole field model -- input and output in SM cartesian coords

    double x_mag[3];
    double b_mag[3], b_car[3];

    // Convert to magnetic dipole coords:
    sm_to_mag_d_(itime_in, x_in, x_mag);
    cardeg(x_mag);

    bmodel_dipole(x_mag, b_mag);

    transform_data_sph2car(x_mag[1], x_mag[2], b_mag, b_car);
    mag_to_sm_d_(itime_in, b_car, b_out);

}


void dipole_geo(int itime_in[2], double x_in[3], double b_out[3]) {
    //  Tilted-dipole magnetic field model.
    //  x_in:   geographic R, lat, lon (R in Earth radii)
    //  b_out:  B-field in nanotesla: Br, Btheta, Bphi
    //          Br: Outward (radial)
    //          Btheta: Southward (colatitude)
    //          Bphi:   Eastward  (longitude)

    double x_in_geocar[3], x_in_magdeg[3];
    double b_magsph[3], b_magcar[3], b_geosph[3], b_geocar[3];

    x_in_geocar = {x_in[0], x_in[1], x_in[2]};

    // Rotate to mag. dipole coords
    degcar(x_in_geocar);
    geo_to_mag_d_(itime_in, x_in_geocar, x_in_magdeg);
    cardeg(x_in_magdeg);

    // Get field in dipole coordinates:
    bmodel_dipole(x_in_magdeg, b_magsph);

    // Map to geographic (spherical):
    transform_data_sph2car(x_in_magdeg[1], x_in_magdeg[2], b_magsph, b_magcar);
    transform_data_mag2geo(itime_in, b_magcar, b_geocar);
    transform_data_car2sph(x_in[1], x_in[2], b_geocar, b_out);    
}




// ---- IGRF MODEL -----------------------------------------------


void init_igrf(int itime_in[2]) {
    // Initialize IGRF for a given time.
    //      recalc_08 sets the various internal parameters for the 
    //      Fortran code, and stores it in a Fortran common block.

    int year = (int)(itime_in[0]/1000);
    int doy  = itime_in[0] - 1000*year;
    int isec = itime_in[1]*1e-3; 
    int hr = (int)(isec/3600);
    int mn = (int)((isec - hr*3600)/60);
    int sc = (int)(isec - hr*3600 - mn*60);
    double vgsex = -400.0; double vgsey = 0; double vgsez = 0; // Default solar wind params

    // Initialize igrf fortran code
    recalc_08_(&year, &doy, &hr, &mn, &sc, &vgsex, &vgsey, &vgsez);
}


void igrf_geo(double x_in[3], double b_out[3]) {
    // x_in:  geographic R, lat, lon  (R in Earth radii)
    // b_out: B-field in nanotesla: Br, Btheta, Bphi
    //          Br: Outward (radial)
    //          Btheta: Southward (colatitude)
    //          Bphi:   Eastward  (longitude)

    double r = x_in[0];
    double theta = D2R*(90. - x_in[1]); // Colatitude in radians
    double phi   = D2R*(x_in[2]);       // East latitude in radians

    igrf_geo_08_(&r, &theta, &phi, b_out, b_out+1, b_out+2);
}

void igrf_mag_cart(int itime_in[2], double x_in[3], double b_out[3], bool recalc) {
    // bool recalc: whether or not to recalculate the IGRF coefficients
    //              (they're constant per itime_in, so don't if it hasn't changed)
    // x_in, b_out: are in magnetic dipole Cartesian coordinates.

    double x_in_geo[3];
    double b_geosph[3], b_geocar[3];

    if (recalc) {
        init_igrf(itime_in);
    }

    // convert mag cartesian to geographic spherical:
    mag_to_geo_d_(itime_in, x_in, x_in_geo);
    cardeg(x_in_geo);    

    // Get IGRF outputs in geographic spherical:
    igrf_geo(x_in_geo, b_geosph);

    // Map to geographic cartesian:
    transform_data_sph2car(x_in_geo[1], x_in_geo[2], b_geosph, b_geocar);

    // Map to magnetic cartesian:
    transform_data_geo2mag(itime_in, b_geocar, b_out);

}



void bmodel(int itime_in[2], double x_in[3], double tsyg_params[10], int model_number, double b_out[3]) {
    // Selectable B-field model in SM Cartesian coordinates.


    double x_gsm[3];
    double b_tmp[3];            // Temporary coords.
    double b_int[3];            // Internal magnetic field (IGRF or dipole)
    // float  b_ext[3];            // external magnetic field (Tsyg solar wind model)
    float beX, beY, beZ;        // External magnetic field (Tsyg solar wind model)
    float tX, tY, tZ;
    float psi;
    // Tsyganenko model params:

    double iopt = 0;               // Dummy input that does absolutely nothing

    // cout << "psi (geopack): " << geopack1_.PSI << "\n";

    // Get internal magnetic field:
    switch(model_number) {
        case 0: 
        // Dipole field
            dipole_sm(itime_in, x_in, b_out);
            // sm_to_gsm_d_(itime_in, b_tmp, b_int);
        break;

        case 1:
        // IGRF model
            init_igrf(itime_in);
            sm_to_gsm_d_(itime_in, x_in, x_gsm);


             // IGRF model uses Tsyganenko's GSW coordinates. These collapse to the standard
            // GSM coordinates if we set the solar wind velocities vgsex = -400, vgsey = vgsez=0.
            // (Default assigned in init_igrf()).
            igrf_gsw_08_(x_gsm, x_gsm + 1, x_gsm + 2, b_int, b_int + 1, b_int + 2);
            gsm_to_sm_d_(itime_in, b_int, b_out);
        break;

        case 2:
        // Tsyganenko model
            init_igrf(itime_in);
            sm_to_gsm_d_(itime_in, x_in, x_gsm);

            psi = float(geopack1_.PSI);

            tX = float(x_gsm[0]); tY = float(x_gsm[1]); tZ = float(x_gsm[2]);
            t04_s_(&iopt, tsyg_params, &psi, &tX, &tY, &tZ, &beX, &beY, &beZ);
            b_int[0] = -1.*double(beX); 
            b_int[1] = -1.*double(beY); 
            b_int[2] = -1.*double(beZ);
            gsm_to_sm_d_(itime_in, b_int, b_out);
        break;
    }
    // Convert nanotesla -> tesla
    b_out[0]*=1e-9;
    b_out[1]*=1e-9;
    b_out[2]*=1e-9;
}


int trace_fieldline(int itime_in[2], double x_in[3], double x_out[TRACER_MAX][3], 
    double b_out[TRACER_MAX], double ds_in, int model_number, double tsyg_params[10]) {
    // x_in: position in SM cartesian coordinates
    // x_out: [nx3] array of field line positions (SM cartesian coordinates)
    // b_out: magnitude of B-field at each position

    double Bo[3] = {0};
    double Bomag;

    double x_tmp[3];
    double x_tmp_mag, x_cur_mag;

    double dx, dy, dz;

    double ds = ds_in;

    double x_cur[3];
    double x_mag[3];    // Magnetic dipole coordinates

    double x_alt;
    x_cur[0] = x_in[0];
    x_cur[1] = x_in[1];
    x_cur[2] = x_in[2]; 
    

    sm_to_mag_d_(itime_in, x_cur, x_mag);
    cardeg(x_mag);

    double Req = x_mag[0]/pow(cos(D2R*x_mag[1]),2);
    double Rdip = 0;

    double dipole_lat_spacing = (2.*fabs(x_mag[1]))/TRACER_MAX;
    if (x_mag[1] < 0) { dipole_lat_spacing*= -1.0;}

    // Determine which direction to step initially:
    //  -- Make one step, and see if we got closer or further away
    bmodel(itime_in, x_cur, tsyg_params, model_number, Bo);

    // Get unit vectors:
    Bomag = norm(Bo, 3);
    dx = Bo[0]/Bomag; dy = Bo[1]/Bomag; dz = Bo[2]/Bomag;

    // Update x_cur:
    x_tmp[0] += dx*ds;
    x_tmp[1] += dy*ds;
    x_tmp[2] += dz*ds;

    ds = (norm(x_tmp, 3) <= norm(x_cur, 3) ? ds_in : -ds_in);
    // cout << "ds is: " << ds << "\n";
    int i=0;
    while (i < TRACER_MAX) {

        // Get B:
        bmodel(itime_in, x_cur, tsyg_params, model_number, Bo);

        // Get unit vectors:
        Bomag = norm(Bo, 3);

        x_out[i][0] = x_cur[0];
        x_out[i][1] = x_cur[1];
        x_out[i][2] = x_cur[2];

        b_out[i] = Bomag;

        if (model_number == 0) {
            // Dipole model; we can calculate this directly:
            sm_to_mag_d_(itime_in, x_cur, x_mag);
            cardeg(x_mag);
            x_mag[1] -= dipole_lat_spacing;
            x_mag[0] = Req*pow(cos(x_mag[1]*D2R),2.);
            degcar(x_mag);
            mag_to_sm_d_(itime_in, x_mag, x_cur);
        } else {

            // Update x_cur:
            dx = Bo[0]/Bomag; dy = Bo[1]/Bomag; dz = Bo[2]/Bomag;

            x_cur[0] += dx*ds;
            x_cur[1] += dy*ds;
            x_cur[2] += dz*ds;
        }
        x_alt = norm(x_cur, 3);
        // cout << x_alt << endl;

        // Analytical dipole calculation:

        Rdip = Req*pow(cos(D2R*x_cur[1]),2);

        // cout << "xcur(" << i << ") : ";
        // sm_to_mag_d_(itime_in, x_cur, x_mag);
        // cardeg(x_mag);
        // print_array(x_mag, 3);


        if (x_alt < 1) {
            // cout << "alt less than 1" << endl;
            break;
        }

        i++;
    }
    // cout << "Stopped after " << i << " steps\n";
    return i;
}

