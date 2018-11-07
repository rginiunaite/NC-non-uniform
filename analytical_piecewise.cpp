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

    double D = 1;//0.05; // to 10^5 \nu m^2/h diffusion coefficient
    double k_reac = 0.105;//0.105; // reaction term
    double C0 = 1.0; // initial chemo concentration in the first part
    double length_x_initial = 1.0; // initial length of the domain
    int solution_grid = 100; // solution grid
    double length_x = length_x_initial; // this will be the actual length
    double beta = 0.5; // initially chemo is non-zero up to beta
    double non_growing_final = beta*0.5; // up to here domain does not grow
    double beta_non_grow = non_growing_final; // initially chemo is non-zero up to beta_non_grow in the non-growing part
    double beta_grow = beta - beta_non_grow; // this is non-zero chemo for the growing part
    double length_constant_initial = (length_x_initial * non_growing_final ); // initial length of constant part
    double length_growing_initial = (length_x_initial - length_x_initial * non_growing_final ); // initial length of growing part
    double length_x_growing_part = (length_x_initial - length_x_initial * non_growing_final); // length of the growing part


    // end of the non-growing domain in grid coordinates

    int gamma = int(non_growing_final*double(solution_grid));

    int t_final = 20; // simulation length
    int n = 1000; // terms fro truncating sum

    MatrixXd an = MatrixXd::Zero(n,solution_grid);
    VectorXd an_temp = VectorXd::Zero(n); // Fourier series coefficients
    MatrixXd bn = MatrixXd::Zero(n,solution_grid); // cos terms
    VectorXd bn_temp = VectorXd::Zero(n);
    MatrixXd cn = MatrixXd::Zero(n,solution_grid); // exponential terms, also dependent on space
    VectorXd cn_temp = VectorXd::Zero(n);
    VectorXd dn = VectorXd::Zero(n);  // all terms multiplied


    MatrixXd chemo = MatrixXd::Zero(solution_grid,t_final+1); // exponential terms
    MatrixXi grid_changes = MatrixXi::Zero(solution_grid,t_final+1); // this will keep track of changing coordinates

    double alpha = 0.1; // exponential growth rate


    VectorXd strain = VectorXd::Zero(solution_grid);


//    // constant
//    for (int i = 0; i < solution_grid; i++) {
//        strain(i) = alpha;// 0;//linear_par * double(theta1) / double(space_grid_controller);//0;
//    }

    // piecewise constant
    for (int i = 0; i < gamma; i++) {
        strain(i) = 0;// 0;//linear_par * double(theta1) / double(space_grid_controller);//0;
    }
    for (int i = gamma; i < solution_grid; i++) {
        strain(i) = alpha;// 0;//linear_par * double(theta1) / double(space_grid_controller);//0;
    }


    // phi function for steady state solution and involved boundary conditions

    VectorXd phi = VectorXd::Zero(solution_grid);
    ///coefficients for phi, which are constant throughout the simulations
    double d1 = D/(length_constant_initial*length_constant_initial);
    double d2 = (D)/(length_constant_initial*length_constant_initial); // this one will change!!!
    double f1 = k_reac;
    double f2 = k_reac - alpha;
    double sqrt_f1d1 = sqrt(f1/d1);
    double sqrt_f2d2 = sqrt(f2/d2);
    double denom1 = - (sqrt_f1d1 * ( exp(-sqrt_f1d1)+ exp(sqrt_f1d1)));
    double denom2 = - (sqrt_f2d2* ( exp(-sqrt_f2d2)+ exp(sqrt_f2d2)));

    // growth function

    // I will mainly use its first derivative with respect to x

    VectorXd Gamma_x = VectorXd::Zero(solution_grid);
    VectorXd Gamma = VectorXd::Zero(solution_grid);

        for (int i = 0; i < gamma; i++) {
            Gamma_x(i) = length_constant_initial * exp(0 * strain(i));

        }
    for (int i = gamma;i < solution_grid; i++) {
        Gamma_x(i) = length_growing_initial * exp(0 * strain(i));
    }



    // an coefficients are independent of space
    double x;

    for (int j = 0; j < gamma; j++) {

        an(0, j) = beta_non_grow * C0 / length_constant_initial;
    }
    //for (int j = int(non_growing_final*double(solution_grid));j < solution_grid; j++) {
    for (int j = 0 ;j < solution_grid- gamma; j++) {
        an(0 , j+ gamma) = beta_grow * C0 / length_growing_initial;
    }

    for (int i = 1; i< n; i++){

        for (int j = 0; j < gamma; j++) {

            an(i,j) = 2*C0 / (i*M_PI) * sin(i *M_PI*beta_non_grow /length_constant_initial);
        }
        for (int j = 0; j < solution_grid - gamma; j++) {
            an(i ,j+ gamma) = 2*C0 / (i*M_PI) * sin(i *M_PI*beta_grow / length_growing_initial);
        }


    }


    for (int t = 0; t <t_final+1;t++){



        // update the strain rate

        length_x_growing_part = length_growing_initial * exp(t * strain(solution_grid-1)); // strain rate is the same on the growing part,so take any coefficient


        //length_x = length_x_initial * exp(alpha*t);

        //cout << "length " << length_x << endl;


        for (int i = 0; i< n;i++){

            for (int xL = 0; xL < gamma; xL++) {
                cn(i,xL) = exp(- (D*i*i*M_PI*M_PI*t)/(length_constant_initial*length_constant_initial) +t*(k_reac) );
                x = double(xL)/solution_grid;
                bn(i,xL) = cos (i*M_PI*x/length_constant_initial);
            }
            for (int xL = 0 ;xL < solution_grid- gamma; xL++) {
                cn(i ,xL+ gamma) = exp(- (D*i*i*M_PI*M_PI*(1.0- exp(-2.0*alpha*t)))/(2.0*alpha*length_constant_initial*length_constant_initial) +t*(k_reac-alpha) );
                x = length_x_growing_part*double(xL)/solution_grid;
                bn(i ,xL + gamma) = cos (i*M_PI*x/length_x_growing_part);
            }
        }



        for (int y = 0; y < solution_grid; y++){

            d2 = (D)/(length_x_growing_part*length_x_growing_part);
            sqrt_f2d2 = sqrt(f1/d1);
            denom2 = - (sqrt_f2d2* ( exp(-sqrt_f2d2)+ exp(sqrt_f2d2)));


            if (t>0) {
                if (y < gamma) {
                    cout << "denom 2 " << denom2 << endl;


                    phi(y) = (chemo(gamma, t - 1) - chemo(gamma - 1, t - 1)) / denom1 *
                             (exp(-sqrt_f1d1 * y) - exp(sqrt_f1d1 * y));
                    //cout << "phy first part " << phi(y) << endl;
                } else {
                    phi(y) = (chemo(gamma - 1, t - 1) - chemo(gamma, t - 1)) / denom1 *
                             (exp(-sqrt_f1d1 * y) - exp(sqrt_f1d1 * y));
                    //cout << "phy second part " << phi(y) << endl;
                }
            }

            // take one by one columns of coefficients corresponding to the a1,a2,a3,...
            an_temp = an.col(y);
            bn_temp = bn.col(y) ;
            cn_temp = cn.col(y);


            //cout << " print the entire column " <<  bn.col(y) << endl;
            //cout << " how it copied " <<  bn_temp << endl;
            // element wise multiplication
            for (int i = 0; i < n;i++){
                dn(i) = an_temp(i)*bn_temp(i)*cn_temp(i);
//                if (y == 0 && t ==0){
//                    cout << "bn before " << bn(i,0) << endl;
//
//                    cout << "an(0) " << an(i) << endl;
//                    cout << "bn(0) " << bn_temp(i) << endl;
//                    cout << "cn(0) " << cn(i) << endl;
//                    cout << "dn(0) " << dn(i) << endl;
//                   cout << "dndn " << dn(i) << endl;
//                }

            }
            chemo(y,t) = dn.sum() + phi(y);

        }



        /*
         * piecewise constant, first constant
         *
         * */

        for (int i = 0; i < int(non_growing_final*double(solution_grid)); i++) {
            grid_changes(i,t) = i; //grid does not changes
        }
        for (int i = 0;  i< solution_grid - gamma; i++){
            grid_changes(i + gamma,t) = i * length_x_growing_part/length_growing_initial + gamma;
        }



    }

    ofstream output("analytical_solution.csv");

    // output << "x, y, z, u" << "\n" << endl;


    for (int i = 0; i < solution_grid; i++) {
        for (int j = 0; j < t_final+1; j++) {
            output << chemo(i, j) << ", ";
        }
        output << "\n" << endl;
    }


    ofstream output2("grid_changes.csv");

    for (int i = 0; i < solution_grid; i++) {
        for (int j = 0; j < t_final+1; j++) {
            output2 << grid_changes(i, j) << ", ";
        }
        output2 << "\n" << endl;
    }

}