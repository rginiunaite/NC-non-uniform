//
// Created by rasa on 27/12/18.
//
/*
 * Reaction-diffusion equation on a non-uniformly growing domain, 1D, explicit method, piece-wise non-uniform growth, piecewise two parts, piecewise three parts with one linear
 * Initially there will be no cell dynamics
 * */



#include "Aboria.h"
#include <random>
#include <map>
#include <Eigen/Core>
#include <algorithm> // std::unique_copy
#include <iostream>// writing on a text file
#include <fstream>
#include <math.h>
#include <assert.h>

using namespace std;
using namespace Aboria;
using namespace Eigen; // objects VectorXd, MatrixXd


int main() {


// parameters

// model parameters

/*
 * Initial length of the domain is 300 \mu m, I will create a matrix of chemoattractant which will be of length 300/10= 30.
 * All the length parameters will be scaled consistently, i.e. /10.
 */


    double space_grid_controller = 1000.0;

    double domain_length = 1.0; //this variable is for the actual domain length, since it will be increasing
    double Lt_old = domain_length;
    int length_x =
            int(domain_length) * int(space_grid_controller); // length in x direction of the chemoattractant matrix
    double initial_domain_length = domain_length;
    const int length_y = 1; // length in y direction of the chemoattractant matrix
    const double final_time = 20.0; // number of timesteps, 1min - 1timestep, from 6h tp 24hours.


// parameters for the dynamics of chemoattractant concentration

    double D = 0.00001;//0.00001;//0.05; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0.0; // initialise time
    double dt = 0.01; // time step
    double dt_init = dt;
    int number_time = int(1 / dt_init); // how many timesteps in 1min, which is the actual simulation timestep
    double dx =
            1.0 / double(space_grid_controller); // space step in x direction, double to be consistent with other types
    double dy = 1; // space step in y direction
    double kai = 0.00001;//0;//0.1 // to 1 /h production rate of chemoattractant


    // reaction rate
    double k_reac = 0;//0.1;//0.105;//0.03;//.205; // reaction term


    // domain growth parameters


    /*
    * strain rate
    * */

    //double linear_par = 0.0001;//05;


    // for comparison with analytical
    double alpha = 0.1;//before 0.1


    VectorXd strain = VectorXd::Zero(length_x);


    // constant
//    for (int i = 0; i < length_x; i++) {
//        strain(i) = alpha;// 0;//linear_par * double(theta1) / double(space_grid_controller);//0;
//    }

    //piecewise constant, two parts

    int theta1 = int(0.5 * space_grid_controller);
    double thetasmall = 0.5;
    //int theta2 = int(0.7 * space_grid_controller);


    // first part it is linear growth
    for (int i = 0; i < theta1; i++) {
        strain(i) = 0;//linear_par * double(theta1) /
        // double(space_grid_controller);//linear_par * (double(i) / double(space_grid_controller));
    }


    // second part is constant
    for (int i = theta1; i < length_x; i++) {
        strain(i) = alpha;// 0.002;//0.5 * linear_par * double(theta1) / double(space_grid_controller); // constant to where it was
        //strain(i,j) = linear_par*theta1/(theta1- (theta2-1))*(i-(theta2-1)); // linearly decreasing, if I do this do not forget to change Gamma
    }




    // three parts, including spatial linear growth

////        int theta1 = int(0.4 * space_grid_controller);
////    int theta2 = int(0.7 * space_grid_controller);
//
    // first part it is linear growth
//    for (int i = 0; i < theta1; i++) {
//        strain(i) = alpha * (double(i) / double(space_grid_controller));
//    }
//
//
//    // second part is constant
//    for (int i = theta1; i < length_x; i++) {
//        strain(i) = alpha * double(theta1) / double(space_grid_controller); // constant to where it was
//        //strain(i,j) = linear_par*theta1/(theta1- (theta2-1))*(i-(theta2-1)); // linearly decreasing, if I do this do not forget to change Gamma
//    }

////    // third part no growth, constant
////    for (int i = theta2; i < length_x; i++) {
////        strain(i) = 0;//linear_par * double(theta1) / double(space_grid_controller);//0;
////    }

    // growth function

    // I will mainly use its first derivative with respect to x

    VectorXd Gamma_x = VectorXd::Zero(length_x);
    VectorXd Gamma = VectorXd::Zero(length_x);
    VectorXd Gamma_t = VectorXd::Zero(length_x);
    VectorXd Gamma_old = VectorXd::Zero(length_x);

    for (int i = 0; i < length_x; i++) {
        Gamma_x(i) = exp(0 * strain(i));
        Gamma(i) = i / space_grid_controller;
        Gamma_old(i) = Gamma(i);
    }

    // for total length
    double Lt = 0;
    double Ltdot = 0;

    /*
    * initialise a matrix that stores values of concentration of chemoattractant
    */

    MatrixXd chemo = MatrixXd::Zero(length_x, length_y);
    MatrixXd chemo_new = MatrixXd::Zero(length_x, length_y);

//    // uniform initial conditions
//    for (int i = 0; i < length_x; i++) {
//        for (int j = 0; j < length_y; j++) {
//            chemo(i, j) = 1; // uniform concentration initially
//            chemo_new(i, j) = 1; // this is for later updates
//        }
//    }e

    // non uniform initial conditions
    double beta = 1.0; // up to here the initial chemo concentration is C_0
    double C0 = 1.0; // initially non-zero, afterwards zero
    double n = 10.0; // for 1 - 0.5 cos(n \pi x)

    for (int i = 0; i < beta * space_grid_controller; i++) {
        for (int j = 0; j < length_y; j++) {
            chemo(i, j) = cos(
                    M_PI * i / space_grid_controller);//1;//C0 - 0.5 * cos( M_PI * i/space_grid_controller * n);
        }
    }

    // initialise internalisation matrix
    MatrixXd intern = MatrixXd::Zero(length_x, length_y);


    // four columns for x, y, z, u (z is necessary for paraview)

    // form a matrix which would store x,y,z,u

    MatrixXd chemo_3col(length_x * length_y, 4), chemo_3col_ind(length_x * length_y,
                                                                2); // need for because that is how paraview accepts data, third dimension is just zeros

    // x, y coord, 1st and 2nd columns respectively
    int k = 0;
    // it has to be 3D for paraview
    while (k < length_x * length_y) {
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                chemo_3col_ind(k, 0) = i;
                chemo_3col_ind(k, 1) = j;
                chemo_3col(k, 2) = 0;
                k += 1;
            }
        }
    }

    // save the x coordinates, scaling only based on the grid
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 0) = chemo_3col_ind(i, 0) / double(space_grid_controller);
    }



    // u column
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
    }


    // y coordinates, 1D so nothing changes
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 1) = chemo_3col_ind(i, 1);
    }


    int counter = 0;


    // save data to plot chemoattractant concentration
    ofstream outputinit("matrix_non_uniform0.csv");

    //output << "x, y, z, u" << "\n" << endl;


    for (int i = 0; i < length_x * length_y; i++) {
        for (int j = 0; j < 4; j++) {
            outputinit << chemo_3col(i, j) << ", ";
        }
        outputinit << "\n" << endl;
    }


    //ofstream output2("track_point" + to_string(t) + ".csv");

    ofstream output2in("track_point0.csv");

    output2in << Gamma(length_x / 2) << endl;

    ofstream output3in("track2_point0.csv");

    output3in << Gamma(int(length_x / 4)) << endl;

    ofstream output4in("track3_point.csv");

    output4in << Gamma(3 * int(length_x / 4)) << endl;



    //for each timestep
    while (t < final_time) {


        t = t + dt;

        counter = counter + 1;

        // update the strain rate
        for (int i = 0; i < length_x; i++) {
            Gamma_x(i) = exp(t * strain(i));
        }


        // domain length


        /*
         * this is important and it will change based on strain rates, now since there is linear growth in the first section,
         * the first factor appears due to integration
         *
         * */
//        if (t != 0) {
//            domain_length = 1.0 / (t * linear_par) * (Gamma_x(theta1 - 1) - 1) +
//                            ((theta2real - 1) - (theta1real - 1)) * Gamma_x(theta2 - 1) + initial_domain_length - 1 -
//                            (theta2real - 1);
//        }


        /*
 * Constant growth, for presentation
 *
 * */


//        for (int i = 0; i < length_x; i++) {
//            for (int j = 0; j < length_y; j++) {
//                Gamma(i) = (double(i) / double(space_grid_controller)) * Gamma_x(i);
//            }
//        }
//



        /*
         * Piecewise constant // all linear, for presentation
         * */


        for (int i = 0; i < theta1; i++) {
            Gamma(i) = (double(i) / double(space_grid_controller)) * (Gamma_x(i)); // linearly increasing

        }



        // second part is constant
        for (int i = theta1; i < length_x; i++) {

            Gamma(i) = Gamma(theta1 - 1) +
                       (double(i) / double(space_grid_controller) - double(theta1 - 1) /
                                                                    double(space_grid_controller)) * Gamma_x(i);
//

            //Gamma(i,j) = ; // linearly decreasing, if I do this do not forget to change Gamma

        }



        /*
 * three different regions
 * */


//        for (int i = 0; i < theta1; i++) {
//            for (int j = 0; j < length_y; j++) {
//                Gamma(i) = 1.0 / (t * alpha) * (Gamma_x(i) - 1); // linearly increasing
//
//            }
//        }
//
//
//        // second part is constant
//        for (int i = theta1; i < length_x; i++) {
//            Gamma(i) = 1.0 / (t * alpha) * (Gamma_x(theta1 - 1) - 1) + (double(i) / double(space_grid_controller) -
//                                                                             double(theta1 - 1) /
//                                                                             double(space_grid_controller)) *
//                                                                            Gamma_x(i); // first was linear, this constant
//            //Gamma(i,j) = ; // linearly decreasing, if I do this do not forget to change Gamma
//
//        }

//
////        // third part no growth, constant
////        for (int i = theta2; i < length_x; i++) {
////            for (int j = 0; j < length_y; j++) {
////                Gamma(i) = 1.0 / (t * alpha) * (Gamma_x(theta1 - 1) - 1) +
////                           (double(theta2 - 1) / double(space_grid_controller) -
////                            double(theta1 - 1) / double(space_grid_controller)) * Gamma_x(theta2 - 1) +
////                           double(i) / double(space_grid_controller) -
////                           (double(theta2 - 1) / double(space_grid_controller)); // linear, constant, zero
////            }
////        }
//





        /*
         * Domain growth
         * */

        Lt = thetasmall + thetasmall * exp(alpha * t);
        Lt = Gamma(length_x - 1);


        //   Ltdot = alpha * thetasmall * exp(alpha * t); // piecewise, exact

        Ltdot = (Lt - Lt_old) / dt; //could be used for both, especially should be used if derivative is unknown
        //    cout << "diff dot " << Ltdot << endl;

        Lt_old = Lt;


        // I need Gamma_t for cos verification as well

        for (int i = 0; i < length_x; ++i) {

            Gamma_t(i) = (Gamma(i) - Gamma_old(i)) / dt;

        }

        Gamma_old = Gamma;




        // coefficient of tridiagonal matrix which is contained in the linear system


        // inner coefficients

        for (int i = 1; i < length_x - 1; ++i) {

            chemo_new(i) = dt * (D * 1.0 / (2.0 * dx * dx * Gamma_x(i)) *
                                 ((1.0 / Gamma_x(i) + 1.0 / Gamma_x(i + 1)) * (chemo(i + 1) - chemo(i)) -
                                  (chemo(i) - chemo(i - 1)) * (1.0 / Gamma_x(i) + 1.0 / Gamma_x(i - 1))) +
                                 chemo(i) * k_reac - strain(i) * chemo(i)) + chemo(i) + dt * (-((Gamma_t(i)) / Lt -
                                                                                                Gamma(i) * Ltdot /
                                                                                                (Lt * Lt)) * M_PI *
                                                                                              sin(M_PI * Gamma(i) /
                                                                                                  Lt) + D * cos(M_PI *
                                                                                                                Gamma(i) /
                                                                                                                Lt) *
                                                                                                        (M_PI / Lt) *
                                                                                                        (M_PI / Lt) -
                                                                                              (k_reac - strain(i)) *
                                                                                              cos(M_PI * Gamma(i) /
                                                                                                  Lt));
        }
        // i = 0, will have to loop over all y when I will move to two dimensions

        chemo_new(0) = dt * (D * 1.0 / (2.0 * dx * dx * Gamma_x(0)) *
                             (((1.0 / Gamma_x(0) + 1.0 / Gamma_x(1)) * (chemo(1) - chemo(0))) -
                              (chemo(0) - chemo(1)) * (1.0 / Gamma_x(0) + 1.0 / Gamma_x(1))) +
                             chemo(0) * (k_reac - strain(0))) + chemo(0) +
                       dt * (-((Gamma_t(0)) / Lt - Gamma(0) * Ltdot / (Lt * Lt)) * M_PI * sin(M_PI * Gamma(0) / Lt) +
                             D * cos(M_PI * Gamma(0) / Lt) *
                             (M_PI / Lt) * (M_PI / Lt) -
                             (k_reac - strain(0)) * cos(M_PI * Gamma(0) / Lt));
       // cout << "chemo_new (0) " << chemo_new(0) << endl;

        chemo_new(length_x - 1) = dt * (D * 1.0 / (2.0 * dx * dx * Gamma_x(length_x - 1)) *
                                        (((1.0 / Gamma_x(length_x - 1) + 1.0 / Gamma_x(length_x - 2)) *
                                          (chemo(length_x - 2) - chemo(length_x - 1))) -
                                         (chemo(length_x - 1) - chemo(length_x - 2)) *
                                         (1.0 / Gamma_x(length_x - 1) + 1.0 / Gamma_x(length_x - 2))) +
                                        chemo(length_x - 1) * k_reac - strain(length_x - 1) * chemo(length_x - 1)) +
                                  chemo(length_x - 1) + dt * (-((Gamma_t(length_x - 1)) / Lt -
                                                                Gamma(length_x - 1) * Ltdot / (Lt * Lt)) * M_PI *
                                                              sin(M_PI * Gamma(length_x - 1) / Lt) +
                                                              D * cos(M_PI * Gamma(length_x - 1) / Lt) *
                                                              (M_PI / Lt) * (M_PI / Lt) -
                                                              (k_reac - strain(length_x - 1)) *
                                                              cos(M_PI * Gamma(length_x - 1) / Lt));


        chemo = chemo_new; // update chemo concentration
        //}




        // save the chemoattractant concentration with properly rescaled coordinates

        int counting_first = 0;
        int counting_final = 0;

        for (int a = 0; a < length_x; a++) {
            counting_first = length_y * a;
            //cout << " Gamma a " << Gamma(a) << endl;
            counting_final = counting_first + length_y;
            for (int k = counting_first; k < counting_final; k++) {
                chemo_3col(k, 0) = Gamma(a);
            }
        }




        //cout << " hello " << endl;


//        // save the chemoattractant concentration with properly rescaled coordinates
//        for (int i = 0; i < length_x * length_y; i++) {
//            chemo_3col(i, 0) = chemo_3col_ind(i, 0) * (domain_length / length_x);
//        }




        // u column
        for (int i = 0; i < length_x * length_y; i++) {
            chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
        }


        // update chemoattractant profile after the domain grew
//        cout << "t " << t << endl;
//        double ine = t;
//        if(ine == t){
//            cout << "is it in here "<< ine << endl;
//
//        }

        if (counter % 100 == 0) {

            //if (t == 1 || t == 10 || t == 20 ) {
            //cout << "heeere " << endl;
            // save data to plot chemoattractant concentration
            //ofstream output("matrix_non_uniform" + to_string(t) + ".csv");

            // output << "x, y, z, u" << "\n" << endl;

            // save data to plot chemoattractant concentration
            ofstream output("matrix_non_uniform" + to_string(t) + ".csv");

            //output << "x, y, z, u" << "\n" << endl;


            for (int i = 0; i < length_x * length_y; i++) {
                for (int j = 0; j < 4; j++) {
                    output << chemo_3col(i, j) << ", ";
                }
                output << "\n" << endl;
            }


            //ofstream output2("track_point" + to_string(t) + ".csv");
//
//            ofstream output2("track_point" + to_string(t) + ".csv");
//
//            output2 << Gamma(length_x / 2) << endl;
//
//            ofstream output3("track2_point" + to_string(t) + ".csv");
//
//            output3 << Gamma(int(length_x / 4)) << endl;
//
//            ofstream output4("track3_point" + to_string(t) + ".csv");
//
//            output4 << Gamma(3 * int(length_x / 4)) << endl;


        }


        //   cout << "t " << t << endl;
//        cout << "domain length " << Gamma_x(length_x-1) << endl;
    }


}