
#include <wipp.h>


// Vector L2 norm.
double l2_norm(vector<double> u) {
    double accum = 0.;
    for (int i = 0; i < u.size(); ++i) {
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


