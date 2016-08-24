#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>


#include <consts.h>
#include <wipp.h>

#include <complex>
#include <cmath>

using namespace std;


double kp_to_pp(double kp){
double pp;

// % pp = kp_to_pp(kp)
// % 
// % Function to convert from kp to PP_L of GCPM model
// % I had to actually look with my eyeballs to calculate these values, but
// % this function interpolates

// % By Daniel Golden (dgolden1 at stanford dot edu) May 2010
// % $Id: kp_to_pp.m 969 2010-05-26 22:04:27Z dgolden $

// Ported to C++ 8.2016 asousa

// X values
const static double kp_vals[] = {3.0, 3.3, 3.7, 4.0, 4.3, 4.7, 5.0, 5.3, 5.7, 6.0, 6.3, 6.7, 7.0, 7.3, 7.7, 8.0};
vector<double> kp_map (kp_vals, kp_vals + sizeof(kp_vals)/sizeof(kp_vals[0]));

// Y values
const static double pp_vals[] = {4.4, 4.3, 4.1, 3.9, 3.8, 3.6, 3.5, 3.4, 3.2, 3.1, 2.9, 2.7, 2.6, 2.4, 2.3, 2.1};
vector<double> pp_map (pp_vals, pp_vals + sizeof(pp_vals)/sizeof(pp_vals[0]));

// Fit coefficients
vector<double> coeffs;

// print_vector(kp_map);
// print_vector(pp_map);

// Do a linear fit on the two vectors
polyfit(kp_map, pp_map, coeffs, 1);
// cout << "coeffs: ";
// print_vector(coeffs);

// Evaluate line @ kp.

return coeffs[0] + coeffs[1]*kp;

}