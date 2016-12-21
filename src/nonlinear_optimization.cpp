#include <wipp.h>

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>


using namespace Eigen;
using namespace std;


double get_d(Matrix<double, 8, 3> P1, Matrix<double, 8, 3> P2, Vector3d targ_point, double nx, double ny, double nz);
double get_S(Matrix<double, 8, 3> P1, Matrix<double, 8, 3> P2, Vector3d targ_point, double nx, double ny, double nz);



// Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
typedef _Scalar Scalar;
enum {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
};
typedef Matrix<Scalar,InputsAtCompileTime,1> InputType;
typedef Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
typedef Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

int m_inputs, m_values;

Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

int inputs() const { return m_inputs; }
int values() const { return m_values; }

};

struct my_functor : Functor<double>
{
my_functor(void): Functor<double>(2,3) {}
Vector3d get_d(VectorXd x);

Matrix<double, 8, 3> P1;
Matrix<double, 8, 3> P2;
Vector3d P3;

double x2;  // The third index value (we're not solving for this one)

int operator()(const VectorXd &x, VectorXd &fvec) const
    // This is the function we want to minimize -- it's the minimum distance
    // between the line between planes P1 and P2 to point P3.
    // N is the interpolation weights, common to both planes at P1 and P2
    {
        double n_x = x[0];
        double n_y = x[1];
        double n_z = x2;

        VectorXd N(8);
        N[0] = (1. - n_x)*(1. - n_y)*(1. - n_z);
        N[1] = n_x*(1. - n_y)*(1. - n_z);
        N[2] = n_y*(1. - n_x)*(1. - n_z);
        N[3] = n_x*n_y*1.*(1. - n_z);
        N[4] = (1. - n_x)*(1. - n_y)*n_z;
        N[5] = n_x*(1. - n_y)*n_z;
        N[6] = n_y*(1. - n_x)*n_z;
        N[7] = n_x*n_y*n_z*1.;

        Vector3d pp1 = N.transpose()*P1;
        Vector3d pp2 = N.transpose()*P2;

        double S = (P3 - pp1).dot(pp2 - pp1)/pow((pp2 - pp1).norm(),2);

        if (S > 1) { S=1;} if (S < 0 ) { S = 0; }
        Vector3d d = pp1 + S*(pp2 - pp1);
        
        Vector3d dvec = d - P3;
        fvec[0] = dvec[0];
        fvec[1] = dvec[1];
        fvec[2] = dvec[2];
        return 0;
    }



};




double get_d(Matrix<double, 8, 3> P1, Matrix<double, 8, 3> P2, Vector3d P3, double nx, double ny, double nz) {
        VectorXd N(8);
        N[0] = (1. - nx)*(1. - ny)*(1. - nz)*1.0;
        N[1] = nx       *(1. - ny)*(1. - nz)*1.0;
        N[2] = (1. - nx)*ny       *(1. - nz)*1.0;
        N[3] = nx       *ny       *(1. - nz)*1.0;
        N[4] = (1. - nx)*(1. - ny)*nz*1.0;
        N[5] = nx       *(1. - ny)*nz*1.0;
        N[6] = (1. - nx)*ny       *nz*1.0;
        N[7] = nx       *ny       *nz*1.0;

        Vector3d pp1 = N.transpose()*P1;
        Vector3d pp2 = N.transpose()*P2;

        double S = (P3 - pp1).dot(pp2 - pp1)/pow((pp2 - pp1).norm(),2);
        Vector3d d = pp1 + S*(pp2 - pp1);
        
        // Vector3d d = get_d(x);

        Vector3d dvec = d - P3;

        return dvec.norm();
}


double get_S(Matrix<double, 8, 3> P1, Matrix<double, 8, 3> P2, Vector3d P3, double nx, double ny, double nz) {
        VectorXd N(8);
        N[0] = (1. - nx)*(1. - ny)*(1. - nz)*1.0;
        N[1] = nx       *(1. - ny)*(1. - nz)*1.0;
        N[2] = (1. - nx)*ny       *(1. - nz)*1.0;
        N[3] = nx       *ny       *(1. - nz)*1.0;
        N[4] = (1. - nx)*(1. - ny)*nz*1.0;
        N[5] = nx       *(1. - ny)*nz*1.0;
        N[6] = (1. - nx)*ny       *nz*1.0;
        N[7] = nx       *ny       *nz*1.0;

        Vector3d pp1 = N.transpose()*P1;
        Vector3d pp2 = N.transpose()*P2;

        double S = (P3 - pp1).dot(pp2 - pp1)/pow((pp2 - pp1).norm(),2);

        return S;
}





void find_crossing(rayT cur_frames[8], rayT prev_frames[8], Vector3d targ_point, double n_f, double *n_x, double *n_y) {
// Solve for the interpolating quantities in latitude and longitude, for which the ray
// passes through the center of the EA segment. PRETTY RAD

// Instantiate the function to be solved:
my_functor functor;

VectorXd x(2);
x(0) = 0.5;
x(1) = 0.5;
// x.setConstant(0.5, 0.5);
// Set parameters:
for (int i=0; i < 8; ++i) {
        functor.P1.row(i)= prev_frames[i].pos;
        functor.P2.row(i)= cur_frames[i].pos;
}

functor.P3 = targ_point;
// cout << functor.P1;
functor.x2 = n_f;

NumericalDiff<my_functor> numDiff(functor);
LevenbergMarquardt<NumericalDiff<my_functor>,double> lm(numDiff);
lm.resetParameters();
lm.parameters.maxfev = 10000;
lm.parameters.ftol = 1.0e-10;
lm.parameters.xtol = 1.0e-10;
int ret = lm.minimize(x);


// HybridNonLinearSolver<my_functor> solver(functor);
// int ret = solver.solveNumericalDiff(x);
// // int ret = solver.hybrd1(x);



// cout << "solved: " << x.transpose() << "\t returned code: " << ret << " ";
// cout << "dist: " << get_d(functor.P1, functor.P2, functor.P3, x.data()[0], x.data()[1], functor.x2) << "\n";

double dd = get_d(functor.P1, functor.P2, functor.P3, x.data()[0], x.data()[1], functor.x2);
double S = get_S(functor.P1, functor.P2, functor.P3, x.data()[0], x.data()[1], functor.x2);
// cout << " dist: " << dd << " S: " << S << "\n";

// Threshold the solution -- In order for the closest point to be between the two
// planes, S must range between 0 and 1.
// Adding a little buffer lets the closest point be outside the planes, but not
// wildly far away. Also confirm that the distance between the closest point and
// the EA segment center is within L_MARGIN (just a good idea to check)
if ( (S >=-2) && ( S <= 2) && (fabs(dd) <= L_MARGIN) ){
    *n_x = x.data()[0];
    *n_y = x.data()[1];
} else {
    *n_x = 20.;
    *n_y = 20.;
}
    


// int main(int argc, char *argv[])
// {
// VectorXd x(2);
// x(0) = 0.5;
// x(1) = 0.5;
// std::cout << "x: " << x << std::endl;

// my_functor functor;

// functor.P1 << 0, 0, 0,
//       1, 0, 0,
//       1, 2, 0,
//       0, 2, 0;

// functor.P2 << 1, 0, 2,
//       2.5,0, 2,
//       2.5, 3, 2,
//       0, 3, 2;
// functor.P3 << 1, 0.2, 1.9;

// NumericalDiff<my_functor> numDiff(functor);

// LevenbergMarquardt<NumericalDiff<my_functor>,double> lm(numDiff);
// lm.parameters.maxfev = 10000;
// lm.parameters.xtol = 1.0e-10;
// std::cout << lm.parameters.maxfev << std::endl;

// int ret = lm.minimize(x);
// // std::cout << lm.iter << std::endl;
// std::cout << ret << std::endl;

// std::cout << "x that minimizes the function: " << x << std::endl;



// return 0;
}