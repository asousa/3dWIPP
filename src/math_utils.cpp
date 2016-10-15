
#include <wipp.h>

using namespace std;
using namespace Eigen;

// Vector L2 norm.
double l2_norm(vector<double> u) {
    double accum = 0.;
    for (int i = 0; i < u.size(); ++i) {
        accum += u[i] * u[i];
    }
    return sqrt(accum);
}

// array L2 norm.
double norm(double u[], int size) {
    double accum = 0.;
    for (int i = 0; i < size; ++i) {
        accum += u[i] * u[i];
    }
    return sqrt(accum);
}


// Multiply each element of a vector u by a scalar v.
vector<double> scalar_multiply(vector<double> u, double v) {
    vector<double> out;
    for (int i = 0; i < u.size(); ++i) {
        out.push_back(u[i]*v);
    }
    return out;
}


// Print out all elements of a vector
void print_vector(vector<double> u) {
    for (int i=0; i < u.size(); i++) {
        cout << u[i] << " ";
    }
    cout << "\n";
}


void print_array(double* arr, double len) {
    for (int i=0; i < len; i++) {
        cout << arr[i] << " ";
    }
    cout << "\n";
}

double dot_product(vector<double>u, vector<double> v) {
    double out =0;
    for (int i =0; i < u.size(); i++) {
        out += u[i]*v[i];
    }
}

vector<double> add(vector<double>u, vector<double> v) {
    vector <double> out;
    out.resize(u.size());
    for (int i =0; i < u.size(); i++) {
        out[i] = u[i] + v[i];
    }
}



// void cardeg(double x[3]) {
// // Cartesian to polar (degrees)
//     carsph(x);
//     x[1] = R2D*x[1];
//     x[2] = R2D*x[2];
// }

// void degcar(double x[3]) {
// // Polar to Cartesian (degrees)
//     x[1] = D2R*x[1];
//     x[2] = D2R*x[2];
//     sphcar(x);
// }


// void carsph(double x[3]) {
//     // in-place rotation from Cartesian to Spherical (radians)
//     // output is R, Theta, Phi

//     double lat, lon, r;
    
//     r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
//     lon = atan2(x[1], x[0]);
//     lat = asin(x[2]/r);

//     x[0] = r;
//     x[1] = lat;
//     x[2] = lon;
// }

// void sphcar(double x[3]) {
//     // in-place rotation from Spherical to Cartesian (radians)
//     // input is R, Theta, Phi, output is x y z
//     double lat, lon, r;
//     r = x[0]; lat = x[1]; lon = x[2];
//     x[0] = r*cos(lat)*cos(lon);
//     x[1] = r*cos(lat)*sin(lon);
//     x[2] = r*sin(lat);
// }

// void transform_data_sphcar(double data[3], double lat, double lon) {
//     // Map a vector field from spherical to cartesian coordinates
//     double M[3][3];
    
//     double d_out[3] = {0};

//     double theta = D2R*(90. - lat);
//     double phi = D2R*lon;

//     double st = sin(theta);
//     double sp = sin(phi);
//     double ct = cos(theta);
//     double cp = cos(phi);

//     // Transformation matrix
//     M[0][0] = st*cp;    M[0][1] = ct*cp;   M[0][2] = -sp;
//     M[1][0] = st*sp;    M[1][1] = ct*sp;   M[1][2] = cp;
//     M[2][0] = ct;       M[2][1] = -st;     M[2][2] = 0;

//     // Matrix multiply
//     for (int row = 0; row < 3; row++) {
//         for (int col = 0; col < 3; col++) {
//             d_out[row] += M[row][col]*data[col];
//         }
//     }

//     // cout << "dout: ";
//     // print_array(d_out, 3);

//     data[0] = d_out[0];
//     data[1] = d_out[1];
//     data[2] = d_out[2];

// }

// void transform_data_carsph(double data[3], double lat, double lon) {
//     // Map a vectyor field from cartesian to spherical coordinates
//  double M[3][3];
    
//     double d_out[3] = {0};

//     double theta = D2R*(90. - lat);
//     double phi = D2R*lon;

//     double st = sin(theta);
//     double sp = sin(phi);
//     double ct = cos(theta);
//     double cp = cos(phi);

//     // Transformation matrix
//     M[0][0] = st*cp;    M[0][1] = st*sp;   M[0][2] = cp;
//     M[1][0] = ct*cp;    M[1][1] = ct*sp;   M[1][2] = -sp;
//     M[2][0] = -st;       M[2][1] = cp;     M[2][2] = 0;

//     // Matrix multiply
//     for (int row = 0; row < 3; row++) {
//         for (int col = 0; col < 3; col++) {
//             d_out[row] += M[row][col]*data[col];
//         }
//     }

//     // cout << "dout: ";
//     // print_array(d_out, 3);

//     data[0] = d_out[0];
//     data[1] = d_out[1];
//     data[2] = d_out[2];

// }


// void transform_data_geomag(int itime_in[2], double d_geo[3], double d_mag[3]) {

//     // Basis vectors in input frame:
//     double A1[3] = {1, 0, 0};
//     double A2[3] = {0, 1, 0};
//     double A3[3] = {0, 0, 1};

//     // Basis vectors in output frame:
//     double B1[3];
//     double B2[3];
//     double B3[3];

//     double data_mag[3];

//     geo_to_mag_d_(itime_in, A1, B1);
//     geo_to_mag_d_(itime_in, A2, B2);
//     geo_to_mag_d_(itime_in, A3, B3);

//     // Inner product: 
//     d_mag[0] = d_geo[0]*B1[0] + d_geo[1]*B2[0] + d_geo[2]*B3[0];
//     d_mag[1] = d_geo[0]*B1[1] + d_geo[1]*B2[1] + d_geo[2]*B3[1];
//     d_mag[2] = d_geo[0]*B1[2] + d_geo[1]*B2[2] + d_geo[2]*B3[2];

// }

