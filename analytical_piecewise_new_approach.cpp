//analytical solution of reaction equation on a piecewise uniformly growing domain. First part up to \beta*0.5 does not grow, the other grows exponentially

#include <numeric>
#include <random>
#include <map>
#include <Eigen/Core>
#include <algorithm> // std::unique_copy
#include <iostream>// writing on a text file
#include <fstream>
#include <math.h>
#include <assert.h>

using namespace std;
using namespace Eigen; // objects VectorXd, MatrixXd
#include <math.h>

int main() {


    // model parameters

    double D = 0.01;//0.05; // to 10^5 \nu m^2/h diffusion coefficient
    double k_reac = 0.105;//0.105; // reaction term
    double C0 = 1.0; // initial chemo concentration in the first part
    double length_x_initial = 1.0; // initial length of the domain
    int solution_grid = 100; // solution grid
    double length_x = length_x_initial; // this will be the actual length
    double beta = 1.0; // initially chemo is non-zero up to beta
    double non_growing_final = beta * 0.5; // up to here domain does not grow
    double beta_non_grow = non_growing_final; // initially chemo is non-zero up to beta_non_grow in the non-growing part
    double beta_grow = beta - beta_non_grow; // this is non-zero chemo for the growing part
    double length_constant_initial = (length_x_initial * non_growing_final); // initial length of constant part
    double length_growing_initial = (length_x_initial -
                                     length_x_initial * non_growing_final); // initial length of growing part
    double length_x_growing_part = (length_x_initial -
                                    length_x_initial * non_growing_final); // length of the growing part
    double length_x_growing_part_derivative = 0.0;

    int t_final = 20; // simulation length
    int n = 1000; // terms fro truncating sum

    MatrixXd an = MatrixXd::Zero(n, solution_grid);
    VectorXd an_temp = VectorXd::Zero(n); // Fourier series coefficients
    MatrixXd bn = MatrixXd::Zero(n, solution_grid); // cos terms
    VectorXd bn_temp = VectorXd::Zero(n);
    MatrixXd cn = MatrixXd::Zero(n, solution_grid); // exponential terms, also dependent on space
    VectorXd cn_temp = VectorXd::Zero(n);
    VectorXd dn = VectorXd::Zero(n);  // all terms multiplied


    MatrixXd chemo = MatrixXd::Zero(solution_grid, t_final + 1); // exponential terms
    MatrixXi grid_changes = MatrixXi::Zero(solution_grid, t_final + 1); // this will keep track of changing coordinates

    double alpha = 0.1; // exponential growth rate


    VectorXd strain = VectorXd::Zero(solution_grid);


//    // constant
//    for (int i = 0; i < solution_grid; i++) {
//        strain(i) = alpha;// 0;//linear_par * double(theta1) / double(space_grid_controller);//0;
//    }

    // piecewise constant
    for (int i = 0; i < int(non_growing_final * double(solution_grid)); i++) {
        strain(i) = 0;// 0;//linear_par * double(theta1) / double(space_grid_controller);//0;
    }
    for (int i = int(non_growing_final * double(solution_grid)); i < solution_grid; i++) {
        strain(i) = alpha;// 0;//linear_par * double(theta1) / double(space_grid_controller);//0;
    }

    // growth function

    // I will mainly use its first derivative with respect to x

    VectorXd Gamma_x = VectorXd::Zero(solution_grid);
    VectorXd Gamma = VectorXd::Zero(solution_grid);

    for (int i = 0; i < int(non_growing_final * double(solution_grid)); i++) {
        Gamma_x(i) = length_constant_initial * exp(0 * strain(i));

    }
    for (int i = int(non_growing_final * double(solution_grid)); i < solution_grid; i++) {
        Gamma_x(i) = length_growing_initial * exp(0 * strain(i));
    }


    // end of the non-growing domain in grid coordinates

    int gamma = int(non_growing_final * double(solution_grid));



    // update the strain rate

    length_x_growing_part = length_growing_initial * exp(1 * strain(solution_grid -
                                                                    1)); // strain rate is the same on the growing part,so take any coefficient

    length_x_growing_part_derivative = length_growing_initial * alpha * exp(1 * strain(solution_grid -
                                                                                       1)); // strain rate is the same on the growing part,so take any coefficient




    double Qn = integral(Q1timescos, 0.0, length_constant_initial, 100);

}
double Q1timescos(double x, double length_x_growing_part_derivative, double D, double length_x_growing_part, double length_constant_initial, double k_reac, int n) {

    return (2 * gamma * length_x_growing_part_derivative +
            2 * D / (length_x_growing_part + length_constant_initial) +
            k_reac * (x * x + 2 * gamma * (length_x_growing_part * length_constant_initial)) *
            cos(n * M_PI * x / length_constant_initial));

};

double Q2timescos(double x,double length_x_growing_part_derivative, double D, double length_x_growing_part, double length_constant_initial, double k_reac, double alpha, int n) {

    return (2 * length_x_growing_part_derivative - 2 * D / (length_x_growing_part + length_constant_initial) +
            (k_reac - alpha) * (-x * x + 2 * (length_x_growing_part + length_constant_initial) *
                                         (length_x_growing_part + length_constant_initial)) *
            cos(n * M_PI * x / length_x_growing_part));

};


double integral(double(*f)(double x), double a, double b, int n) {
    double step = (b - a) / n; //
    double area = 0.0; // signed area
    for (int i = 0; i < n; i++) {
        area += f(a + (i + 0.5) * step) * step;
    }
};
