#include <wipp.h>


using namespace std;
using namespace Eigen;


// ---- TILTED DIPOLE MODEL --------------------------------------


void bmodel_dipole(double* x_in, double* B_out) {
    // x_in: R, lat, lon in magnetic dipole coordinates (spherical)
    // B_out: Br, Blat, Blon in magnetic dipole coordinates (Upward, Southward, Eastward)
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
