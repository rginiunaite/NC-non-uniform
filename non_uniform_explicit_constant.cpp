//
// Created by Rasa on 18.10.8.
//
/*
 * Reaction-diffusion equation on a non-uniformly growing domain, 2D, explicit method, actually it is uniform growth in this case
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
 *
 */
    double space_grid_controller = 100;

    double domain_length = 1.0; //this variable is for the actual domain length, since it will be increasing
    int length_x = int(domain_length) * int(space_grid_controller); // length in x direction of the chemoattractant matrix
    double initial_domain_length = domain_length;
    const int length_y = 1; // length in y direction of the chemoattractant matrix
    const double final_time = 20.0; // number of timesteps, 1min - 1timestep, from 6h tp 24hours.


// parameters for the dynamics of chemoattractant concentration

    double D = 0.01;//0.05; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0.0; // initialise time
    double dt = 1; // time step
    double dt_init = dt;
    int number_time = int(1 / dt_init); // how many timesteps in 1min, which is the actual simulation timestep
    double dx = 1.0 / double(space_grid_controller); // space step in x direction, double to be consistent with other types
    double dy = 1; // space step in y direction
    double kai = 0.00001;//0;//0.1 // to 1 /h production rate of chemoattractant


    // reaction rate
    double k_reac = 0.2;//.205; // reaction term






    // domain growth parameters

    // for comparison with analytical
    double alpha = 0.1;



    VectorXd strain = VectorXd::Zero(length_x);


    // third part no growth, constant
    for (int i = 0; i < length_x; i++) {
        strain(i) = alpha;// 0;//linear_par * double(theta1) / double(space_grid_controller);//0;
    }

    // growth function

    // I will mainly use its first derivative with respect to x

    VectorXd Gamma_x = VectorXd::Zero(length_x);
    VectorXd Gamma = VectorXd::Zero(length_x);

    for (int i = 0; i < length_x; i++) {
        Gamma_x(i) = exp(0 * strain(i));
    }


    /*
    * initialise a matrix that stores values of concentration of chemoattractant
    */

    MatrixXd chemo = MatrixXd::Zero(length_x, length_y);
    MatrixXd chemo_new = MatrixXd::Zero(length_x, length_y);

    // uniform initial
    for (int i = 0; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            chemo(i, j) = 1; // uniform concentration initially
            chemo_new(i, j) = 1; // this is for later updates
        }
    }


//    // non uniform initial conditions
//    double beta = 0.2; // up to here the initial chemo concentration is C_0
//    double C0 = 0.5; // initially non-zero, afterwards zero
//
//    for (int i = 0; i < beta*space_grid_controller; i++) {
//        for (int j = 0; j < length_y; j++) {
//            chemo(i, j) = C0;
//        }
//    }



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


    // save the chemoattractant concentration with properly rescaled coordinates
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 0) = chemo_3col_ind(i, 0) / double(space_grid_controller);
    }


    // u column
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
    }

    // y and x (initially) column
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 1) = chemo_3col_ind(i, 1);
    }

    // save data to plot chemoattractant concentration
    ofstream output("explicit_matrix_non_uniform" + to_string(0) + ".csv");

    //output << "x, y, z, u" << "\n" << endl;


    for (int i = 0; i < length_x * length_y; i++) {
        for (int j = 0; j < 4; j++) {
            output << chemo_3col(i, j) << ", ";
        }
        output << "\n" << endl; }


    int counter = 0;


    //for each timestep
    while(t< final_time){

        t =t + dt;

        counter = counter + 1;

        // update the strain rate
        for (int i = 0; i < length_x; i++) {
            Gamma_x(i) = exp(t * strain(i));
        }


        /*
         * Dynamics of chemoattractant on a growing domain
         * */



        // inner part, explicit method




        if (length_y >1){
            for (int i = 1; i < length_x - 1; i++) {
                for (int j = 1; j < length_y-1; j++) {
                    chemo_new(i, j) = dt * (D * 1.0 / (Gamma_x(i) *dx*2.0) *
                                            (((1.0 / Gamma_x(i) + 1.0 / Gamma_x(i + 1))) *
                                             (chemo(i + 1, j) - chemo(i, j)) / dx -
                                             (1.0 / Gamma_x(i) + 1.0 / Gamma_x(i - 1)) *
                                             (chemo(i, j) - chemo(i - 1, j)) / dx) + D * (chemo(i,j+1) - 2 * chemo(i,j) + chemo(i,j-1))/dy*dy - strain(i)*chemo(i,j))+ chemo(i,j);
                }

            }
        }

        if (length_y == 1){
            for (int i = 1; i < length_x - 1; i++) {
                chemo_new(i, 0) = dt * (D * 1.0 / (Gamma_x(i) *dx*2.0) *
                                        (((1.0 / Gamma_x(i) + 1.0 / Gamma_x(i + 1))) *
                                         (chemo(i + 1, 0) - chemo(i, 0)) / dx -
                                         (1.0 / Gamma_x(i) + 1.0 / Gamma_x(i - 1)) *
                                         (chemo(i, 0) - chemo(i - 1, 0)) / dx)- strain(i)*chemo(i,0))+ chemo(i,0);
            }
        }

        //boundary, zero flux explicit method


        if (length_y>1) { //if more than 1D
            for (int i = 0; i < length_x; i++) {
                chemo_new(i, 0) = chemo_new(i, 1);
                chemo_new(i, length_y - 1) = chemo_new(i, length_y - 2);
            }

            for (int j = 0; j < length_y; j++) {
                chemo_new(0, j) = chemo_new(1, j);
                chemo_new(length_x - 1, j) = chemo_new(length_x - 1, j);
            }
        }
        else{// if 1D
            chemo_new(0,0) = chemo_new(1,0);
            chemo_new(length_x-1,0) = chemo_new(length_x-2,0);
        }




        chemo = chemo_new; // update chemo concentration
        //}




        /*
         * Constant growth
         *
         * */


        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                Gamma(i) = (double(i) / double(space_grid_controller)) * Gamma_x(i);
            }
        }




        // save the chemoattractant concentration with properly rescaled coordinates

        int counting_first = 0;
        int counting_final = 0;

        for (int a = 0; a < length_x; a++){
            counting_first = length_y*a;
            //cout << " Gamma a " << Gamma(a) << endl;
            counting_final = counting_first + length_y;
            for (int k = counting_first; k < counting_final;k++){
                chemo_3col(k,0) = Gamma(a);
            }
        }


        // u column
        for (int i = 0; i < length_x * length_y; i++) {
            chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
        }


        // update chemoattractant profile after the domain grew

        //if (counter % 5 == 0) {

        //if (t == 400 || t == 800 || t == 1200 || t == 1600) {
            //cout << "heeere " << endl;
            // save data to plot chemoattractant concentration
            ofstream output("explicit_matrix_non_uniform" + to_string(t) + ".csv");

            // output << "x, y, z, u" << "\n" << endl;


            for (int i = 0; i < length_x * length_y; i++) {
                for (int j = 0; j < 4; j++) {
                    output << chemo_3col(i, j) << ", ";
                }
                output << "\n" << endl;
            }
        //}
    }


}