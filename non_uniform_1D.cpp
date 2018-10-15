//
// Created by Rasa on 18.10.8.
//
/*
 * Reaction-diffusion equation on a non-uniformly growing domain, 1D, implicit method
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


int main(){


// parameters

// model parameters

/*
 * Initial length of the domain is 300 \mu m, I will create a matrix of chemoattractant which will be of length 300/10= 30.
 * All the length parameters will be scaled consistently, i.e. /10.
 *
 */


    int space_grid_controller = 1;

    int length_x = (30+1)* space_grid_controller;; // length in x direction of the chemoattractant matrix
    double domain_length = 30+1; //this variable is for the actual domain length, since it will be increasing
    double initial_domain_length = domain_length;
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

    double D = 1; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0; // initialise time
    double dt = 1; // time step
    double dt_init = dt;
    int number_time = int(1/dt_init); // how many timesteps in 1min, which is the actual simulation timestep
    double dx = space_grid_controller; // space step in x direction, double to be consistent with other types
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

    int theta1real = 10;
    int theta2real = 20;
    int theta1 = theta1real * space_grid_controller; // this is for numerical simulations
    int theta2 = theta2real*space_grid_controller;




    double linear_par = 0.001;//05;

    double domain_len_der = 0; // initialise derivative of the domain growth function

    VectorXf strain = VectorXf::Zero(length_x);

    // first part it is linear growth
    for (int i = 0; i < theta1; i++) {
            strain(i) = linear_par * (i);

    }


    // second part is constant
    for (int i = theta1; i < theta2; i++) {
            strain(i) = linear_par*theta1; // constant to where it was
            //strain(i,j) = linear_par*theta1/(theta1- (theta2-1))*(i-(theta2-1)); // linearly decreasing, if I do this do not forget to change Gamma
    }

    // third part no growth, constant
    for (int i = theta2; i < length_x; i++) {
            strain(i) = 0;
    }

    // growth function

    // I will mainly use its first derivative with respect to x

    VectorXf Gamma_x = VectorXf::Zero(length_x);
    VectorXf Gamma = VectorXf::Zero(length_x);

    for (int i = 0; i < length_x; i++) {
            Gamma_x(i) =  exp(0 * strain(i));
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

    // save the chemoattractant concentration, no scaling initially
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 0) = chemo_3col_ind(i, 0);
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
    ofstream output("matrix_non_uniform" + to_string(0) + ".csv");

    //output << "x, y, z, u" << "\n" << endl;


    for (int i = 0; i < length_x * length_y; i++) {
        for (int j = 0; j < 4; j++) {
            output << chemo_3col(i, j) << ", ";
        }
        output << "\n" << endl; }


    //for each timestep
    while (t < final_time){

        t=t + dt;


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
        if (t != 0 ){
            domain_length = 1.0/(t*linear_par)* (Gamma_x(theta1real-1) - 1 ) + ((theta2real-1) - (theta1real-1)) * Gamma_x(theta2-1) + initial_domain_length-1-(theta2real-1);
        }



        /*
         * Domain growth
         * */

        // coefficient of tridiagonal matrix which is contained in the linear system

        // only works if y = 1
        MatrixXf ai = MatrixXf::Zero(length_x,length_y);
        MatrixXf bi = MatrixXf::Zero(length_x,length_y);
        MatrixXf gi = MatrixXf::Zero(length_x,length_y);
        MatrixXf di = MatrixXf::Zero(length_x,length_y);


        // inner coefficients

        for (int i = 1; i < length_x-1; i++){
            for (int j = 0; j< length_y; j++ ){
                ai(i,j) = - D * dt/(2.0*Gamma_x(i)*dx*dx)* (1.0/Gamma_x(i) + 1.0/Gamma_x(i-1));

            }
        }


        for (int i = 1; i < length_x-1; i++){
            for (int j = 0; j< length_y; j++ ){
                bi(i,j) =( 1 + dt*strain(i) + D * dt/(2.0*Gamma_x(i)*dx*dx)*( 1.0/ Gamma_x(i) + 1.0/Gamma_x(i+1)) + D * dt/(2.0*Gamma_x(i)*dx*dx)*(1.0/Gamma_x(i) + 1.0/Gamma_x(i-1))) ;
            }
        }


        for (int i = 1; i < length_x-1; i++){
            for (int j = 0; j< length_y; j++ ){
                gi(i,j) = - D * dt/(2.0*Gamma_x(i)*dx*dx)* (1.0/Gamma_x(i) + 1.0/Gamma_x(i+1));

            }
        }


        // zero flux boundary
        for (int j=0;j<length_y;j++){
            bi(0,j) = 1;
            gi(0,j) = -1;
            ai(length_x-1,j) = -1;
            bi(length_x-1, j) = 1;
        }

        // RHS of linear system Ax = d
        di(0,0) = 0;
        di(length_x - 1, 0) = 0;

        for (int i=1;i<length_x-1;i++){
            di(i,0) = chemo(i,0);
        }


        // for Thomas algorithm I reformulate the linear matrix by doing appropriate changes

        // new matrix coefficients

        MatrixXf cidash = MatrixXf::Zero(length_x,length_y);
        MatrixXf didash = MatrixXf::Zero(length_x,length_y);

        // first entry
        cidash(0,0) = gi(0,0)/bi(0,0);

        for (int i = 1; i < length_x-1; i++){
            for (int j = 0; j< length_y; j++ ){
                cidash(i,j) = gi(i,j)/(bi(i,j)-ai(i,j)*cidash(i-1,j));
            }
        }

        didash(0,0) = di(0,0)/bi(0,0);

        for (int i = 1; i < length_x; i++){
            for (int j = 0; j< length_y; j++ ){

                didash(i,j) = (di(i,j) - ai(i,j)*didash(i-1,j))/(bi(i,j)-ai(i,j)*cidash(i-1,j));


            }
        }


        // Using this decomposition we get that

        chemo_new(length_x-1,0) = didash(length_x-1,0); // since only one entry in y, would need to change if I had two dimensions

        for (int k = length_x - 2; k >= 0; k--){
            for(int j=0; j< length_y; j++){
                chemo_new(k,j) = didash(k,j) - cidash(k,j) * chemo_new (k+1,j);
            }
        }



        chemo = chemo_new; // update chemo concentration
        //}


        for (int i = 0; i < theta1; i++) {
            for (int j = 0; j < length_y; j++) {
                Gamma(i) = 1.0/(t*linear_par) * (Gamma_x(i) - 1);

            }
        }


        // second part is constant
        for (int i = theta1; i < theta2; i++) {
                Gamma(i) = 1.0/(t*linear_par) * (Gamma_x(theta1-1)-1) + (i - (theta1real-1)) * Gamma_x(i); // constant to where it was
                //Gamma(i,j) = ; // linearly decreasing, if I do this do not forget to change Gamma


        }


        // third part no growth, constant
        for (int i = theta2; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                Gamma(i) = 1.0/(t*linear_par) * (Gamma_x(theta1-1) - 1) + (theta2real-1 - (theta1real-1)) * Gamma_x(theta2 -1) + i - (theta2real-1);
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



        // save data to plot chemoattractant concentration
        ofstream output("matrix_non_uniform" + to_string(t) + ".csv");

        // output << "x, y, z, u" << "\n" << endl;


        for (int i = 0; i < length_x * length_y; i++) {
            for (int j = 0; j < 4; j++) {
                output << chemo_3col(i, j) << ", ";
            }
            output << "\n" << endl; }





    }


}