
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

int nearest(double arr[], int arr_len, double target, bool reverse_order) {
    // Index of element in array closest to target. 
    int ind;

    // Find the first value over
    if (reverse_order) {
        ind = distance(arr, lower_bound(arr, arr + arr_len, target, descending_order));
    } else {
        ind = distance(arr, lower_bound(arr, arr + arr_len, target));
    }

    // Is the last value under closer?
    if ( abs(arr[ind - 1] - target) < abs(arr[ind] - target)) { ind--; }  

    // Are we at the end? 
    if (ind == arr_len) { ind--; }

    return ind;
}


double signof(double val) {
    return (0 < val) - (val < 0);
}


bool descending_order(double a, double b) {
// Crushingly simple method to search with lower_bound in descending order
    return a >= b; 
}



// double area3D_polygon(int n, double V[][3]) {

//     Vector3d v1;
//     Vector3d v2;

//     Vector3d cross = {0};

//     for (int x=0; x < n; ++x) {
//         v1 = Map<Vector3d>(V[ x ],       3,1);
//         v2 = Map<Vector3d>(V[mod(x+1,n)],3,1);
//         cross += v1.cross(v2);
//     }

//     return cross.norm()/2.;
// }

