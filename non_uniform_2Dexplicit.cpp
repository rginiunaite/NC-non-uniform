//
// Created by Rasa on 18.10.8.
//
/*
 * Reaction-diffusion equation on a non-uniformly growing domain, 2D, explicit method
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
using namespace Eigen; // objects VectorXf, MatrixXf


int main() {


// parameters

// model parameters

/*
 * Initial length of the domain is 300 \mu m, I will create a matrix of chemoattractant which will be of length 300/10= 30.
 * All the length parameters will be scaled consistently, i.e. /10.
 *
 */


    int length_x = 30; // length in x direction of the chemoattractant matrix
    double domain_length = 30; //this variable is for the actual domain length, since it will be increasing
    double old_length = 30;// this is important for the update of the positions of cells
    const int length_y = 1; // length in y direction of the chemoattractant matrix
    double cell_radius = 0.75;//0.5; // radius of a cell
    const double diameter =
            2 * cell_radius; // diameter of a cell
    const int final_time = 200; // number of timesteps, 1min - 1timestep, from 6h tp 24hours.
    const size_t N = 5; // initial number of cells
    double l_filo_y = 2.75;//2; // sensing radius, filopodia + cell radius
    double l_filo_x = 2.75; // sensing radius, it will have to be rescaled when domain grows
    double l_filo_x_in = l_filo_x; // this value is used for rescaling when domain grows based on initial value
    double l_filo_max = 4.5; // this is the length when two cells which were previously in a chain become dettached
//double diff_conc = 0.1; // sensing threshold, i.e. how much concentration has to be bigger, so that the cell moves in that direction
    int freq_growth = 1; // determines how frequently domain grows (actually not relevant because it will go every timestep)
    int insertion_freq = 1; // determines how frequently new cells are inserted, regulates the density of population
    double speed_l = 0.08;// 0.05;//1;//0.05; // speed of a leader cell
    double increase_fol_speed = 1.3;
    double speed_f = increase_fol_speed * speed_l;//0.05;//0.1;//0.08; // speed of a follower cell
    double chemo_leader = 0.9; //0.5; // phenotypic switching happens when the concentration of chemoattractant is higher than this (presentation video 0.95), no phenotypic switching
    double eps = 1; // for phenotypic switching, the distance has to be that much higher
    const int filo_number = 1; // number of filopodia sent
    int same_dir = 0; // number of steps in the same direction +1, because if 0, then only one step in the same direction
    bool random_pers = true; // persistent movement also when the cell moves randomly
    int count_dir = 0; // this is to count the number of times the cell moved the same direction, up to same_dir for each cell


// distance to the track parameters
    double dist_thres = 1;
    int closest_time;
    int leader_track;





// parameters for the dynamics of chemoattractant concentration

    double D = 0.00001;//0.05; // to 10^5 \nu m^2/h diffusion coefficient /// FOR NO DIFFUSION VIDEO I HAD 0.00001
    double t = 0; // initialise time, redundant
    double dt = 0.5; // time step
    double dt_init = dt;
    int number_time = int(1 / dt_init); // how many timesteps in 1min, which is the actual simulation timestep
    double dx = 1; // space step in x direction, double to be consistent with other types
    double dy = 1; // space step in y direction
    double kai = 0.00001;//0;//0.1 // to 1 /h production rate of chemoattractant


// parameters for internalisation

    double R = cell_radius;//7.5/10; // \nu m cell radius
    double lam = 0.00035;//(100)/10; // to 1000 /h chemoattractant internalisation



    // domain growth parameters


    /*
    * strain rate
    * */

    // regions where the domain growth rate is different
    int theta1 = 10;
    int theta2 = 20;

    double linear_par = 0.001;

    double domain_len_der = 0; // initialise derivative of the domain growth function

    MatrixXf strain = MatrixXf::Zero(length_x, length_y);

    // first part it is linear growth
    for (int i = 0; i < theta1; i++) {
        for (int j = 0; j < length_y; j++) {
            strain(i, j) = linear_par * (i + 1);
        }
    }


    // second part is constant
    for (int i = theta1; i < theta2; i++) {
        for (int j = 0; j < length_y; j++) {
            strain(i, j) = linear_par * theta1;// constant, continuous with the previous
            //strain(i, j) = linear_par * theta1 / (theta1 - (theta2 - 1)) * (i - (theta2 - 1)); // linearly decreasing
        }
    }

    // third part no growth, constant
    for (int i = theta2; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            strain(i, j) = 0;
        }
    }

    // growth function

    // I will mainly use its first derivative with respect to x

    MatrixXf Gamma_x = MatrixXf::Zero(length_x, length_y);

    for (int i = 0; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            Gamma_x(i, j) = i * exp(0 * strain(i, j));
        }
    }


    /*
    * initialise a matrix that stores values of concentration of chemoattractant
    */

    MatrixXf chemo = MatrixXf::Zero(length_x, length_y);
    MatrixXf chemo_new = MatrixXf::Zero(length_x, length_y);

    for (int i = 0; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            chemo(i, j) = 1; // uniform concentration initially
            chemo_new(i, j) = 1; // this is for later updates
        }
    }

    // initialise internalisation matrix
    MatrixXf intern = MatrixXf::Zero(length_x, length_y);


    // four columns for x, y, z, u (z is necessary for paraview)

    // form a matrix which would store x,y,z,u

    MatrixXf chemo_3col(length_x * length_y, 4), chemo_3col_ind(length_x * length_y,
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
        chemo_3col(i, 0) = chemo_3col_ind(i, 0) * (domain_length / length_x);
    }


    // u column
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
    }

    // y and x (initially) column
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 1) = chemo_3col_ind(i, 1);
        chemo_3col(i, 0) = chemo_3col_ind(i, 0);
    }

    // save data to plot chemoattractant concentration
    ofstream output("matrix_non_uniform2D" + to_string(0) + ".csv");

    //output << "x, y, z, u" << "\n" << endl;


    for (int i = 0; i < length_x * length_y; i++) {
        for (int j = 0; j < 4; j++) {
            output << chemo_3col(i, j) << ", ";
        }
        output << "\n" << endl; }


    //for each timestep
    while(t< final_time){

        t =t + dt;

        // update the strain rate
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                cout << "strain " << strain(i, j) << endl;
                Gamma_x(i, j) = exp(t * strain(i, j));
                cout << "Gamma_x " << Gamma_x(i, j) << endl;
            }
        }


        // domain length

        cout << "Gamma last" << Gamma_x(length_x - 1, 0) << endl;

        /*
         * this is important and it will change based on strain rates, now since there is linear growth in the first section,
         * the first factor appears due to integration
         *
         * */
        if (t != 0) {
            domain_length = 1.0 / (t * linear_par) * (Gamma_x(theta1 - 1, 0) - 1) +
                            (theta2 - 1 - (theta1 - 1)) * Gamma_x(theta2 - 1, 0) + length_x - 1 - (theta2 - 1);
        }

        cout << "domain_length" << domain_length << endl;


        /*
         * Dynamics of chemoattractant on a growing domain
         * */



        // inner part, explicit method

        for (int i = 1; i < length_x - 1; i++) {
            for (int j = 0; j < length_y; j++) {
                chemo_new(i, j) = dt * (D * 1.0 / Gamma_x(i, j) / (dx*2.0) *
                                    (((1.0 / Gamma_x(i, j) + 1.0 / Gamma_x(i + 1, j))) *
                                     (chemo(i + 1, j) - chemo(i, j)) / dx -
                                     (1.0 / Gamma_x(i, j) + 1.0 / Gamma_x(i - 1, j)) *
                                    (chemo(i, j) - chemo(i - 1, j)) / dx) - strain(i,j)*chemo(i,j))+ chemo(i,j);
            }

        }

        //boundary, zero flux implicit method

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



        // save the chemoattractant concentration with properly rescaled coordinates
        for (int i = 0; i < length_x * length_y; i++) {
            chemo_3col(i, 0) = chemo_3col_ind(i, 0) * (domain_length / length_x);
        }


        // u column
        for (int i = 0; i < length_x * length_y; i++) {
            chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
        }


        // update chemoattractant profile after the domain grew





        // save data to plot chemoattractant concentration
        ofstream output("matrix_non_uniform2D" + to_string(t) + ".csv");

        //output << "x, y, z, u" << "\n" << endl;


        for (int i = 0; i < length_x * length_y; i++) {
            for (int j = 0; j < 4; j++) {
                output << chemo_3col(i, j) << ", ";
            }
            output << "\n" << endl;
        }


    }


}