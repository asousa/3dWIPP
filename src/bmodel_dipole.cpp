#include <wipp.h>


using namespace std;
using namespace Eigen;

void bmodel_dipole(double* x_in, double* B_out) {
    // x_in: R, lat, lon in magnetic dipole coordinates
    // B_out: Br, Blat, Blon in magnetic dipole coordinates
    const double Bo = 3.12e-5;

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
    B_out[1] = Btheta;      // North 
    B_out[2] = Bphi;        // East

}


int trace_fieldline(double x_in[3], double x_out[TRACER_MAX][3], double ds_in) {
// x_in: R, Lat, Lon in magnetic dipole coordinates (earth radii)
// x_out: [nx3] array of field line coordinates

    double Bo[3] = {0};
    double Bomag, xcurmag;

    double dx, dy, dz;

    double ds = ds_in;


    double x_cur[3];

    x_cur[0] = x_in[0];
    x_cur[1] = x_in[1];
    x_cur[2] = x_in[2]; 
    

    // Determine which direction to step initially:
    // Loosely, negative direction in northern hemisphere.
    // (This might be too-simplistic near the equator with complex field models)
    ds = (x_in[1] > 0 ? -ds_in : ds);

    int i=0;
    while (i < TRACER_MAX) {
    // for (int i=0; i < TRACER_MAX; i++ ) {

        bmodel_dipole(x_cur, Bo);
        transform_data_sphcar(Bo, x_cur[1], x_cur[2]);

        Bomag = sqrt(Bo[0]*Bo[0] + Bo[1]*Bo[1] + Bo[2]*Bo[2]);
        dx = Bo[0]/Bomag; dy = Bo[1]/Bomag; dz = Bo[2]/Bomag;
        
        degcar(x_cur);

        x_cur[0] += dx*ds;
        x_cur[1] += dy*ds;
        x_cur[2] += dz*ds;

        cardeg(x_cur);

        x_out[i][0] = x_cur[0];
        x_out[i][1] = x_cur[1];
        x_out[i][2] = x_cur[2];

        // cout << "xcur(" << i << ") : ";
        // print_array(x_cur,3);

        if (x_cur[0] < 1) {
            cout << "Stopped after " << i << " steps\n";
            break;
        }

        i++;
    }

    return i;

}