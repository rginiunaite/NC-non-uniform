//analytical solution of reaction equation on a uniformly growing domain

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
    double k_reac = 0.105; // reaction term
    double C0 = 1.0; // initial chemo concentration in the first part
    double length_x_initial = 1.0; // initial length of the domain
    int solution_grid = 100; // solution grid
    double length_x = length_x_initial; // this will be the actual length
    double alpha = 0.1; // exponential growth rate
    double beta = 0.2; // initially chemo is non-zero up to beta

    int t_final = 20; // simulation length
    int n = 1000; // terms fro truncating sum

    VectorXd an = VectorXd::Zero(n); // Fourier series coefficients
    MatrixXd bn = MatrixXd::Zero(n,solution_grid); // cos terms
    VectorXd bn_temp = VectorXd::Zero(n);
    VectorXd cn = VectorXd::Zero(n); // exponential terms
    VectorXd dn = VectorXd::Zero(n);  // all terms multiplied


    MatrixXd chemo = MatrixXd::Zero(solution_grid,t_final+1); // exponential terms
    MatrixXi grid_changes = MatrixXi::Zero(solution_grid,t_final+1); // this will keep track of changing coordinates

    double x;

    an(0) = beta*C0 / length_x_initial;

    for (int i = 1; i< n; i++){

        an(i) = 2*C0 / (i*M_PI) * sin(i *M_PI*beta / length_x_initial);
       // cout << "an " << an(i) << endl;
    }


    for (int t = 0; t <t_final+1;t++){

        length_x = length_x_initial * exp(alpha*t);

        //cout << "length " << length_x << endl;

        for (int i = 0; i< n;i++){

            cn(i) = exp(- (D*i*i*M_PI*M_PI*(1.0- exp(-2.0*alpha*t)))/(2.0*alpha*length_x_initial*length_x_initial) +t*(k_reac-alpha) );
            for (int xL= 0;xL< solution_grid;xL++){

                x = length_x*double(xL)/solution_grid;
                bn(i,xL) = cos (i*M_PI*x/length_x);

            }
        }



        for (int y = 0; y < solution_grid; y++){


            bn_temp =   bn.col(y) ;


            //cout << " print the entire column " <<  bn.col(y) << endl;
            //cout << " how it copied " <<  bn_temp << endl;
            // element wise multiplication
            for (int i = 0; i < n;i++){
                dn(i) = an(i)*bn_temp(i)*cn(i);
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
            chemo(y,t) = dn.sum();

        }

        for (int i = 0;i < solution_grid;i++){
            grid_changes(i,t) = i * length_x/length_x_initial;
        }
//        cout << "t " << t << endl;
//        cout << "length_x " << length_x << endl;
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