//
// Created by giniunaite on 07/06/19.
// For domain growth
//


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




VectorXi proportions(double diff_conc, int n_seed) {
//int  main(){

// parameters


    //specify whether first or final part of the domain grows faster

    bool first_part_grows = true; //false if final part grows faster


    double space_grid_controller = 100.0;

    double domain_length = 3.42; //this variable is for the actual domain length, since it will be increasing
    double Lt_old = domain_length;
    int length_x =
            int(domain_length * space_grid_controller); // length in x direction of the chemoattractant matrix
    double initial_domain_length = domain_length;
    const int length_y = int(1.2 * space_grid_controller); // length in y direction of the chemoattractant matrix
    const double final_time = 54; // number of timesteps, 1min - 0.05, now dt =0.01, for 18hours we have 54.
    double final_length = 1014;//real 1014

// parameters for the dynamics of chemoattractant concentration

    double D = 2.0;//0.00001;//0.05; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0.0; // initialise time
    double dt = 0.01; // time step
    double dt_init = dt;
    int number_time = int(1 / dt_init); // how many timesteps in 1min, which is the actual simulation timestep
    double dx = 1.0;// / space_grid_controller; // space step in x direction, double to be consistent with other types
    double dy = 1.0;// / space_grid_controller; // space step in y direction

    // reaction rate
    double k_reac = 1.0; // reaction term


    // cell parameters

    double cell_radius = 7.5;//0.5; // radius of a cell
    const double diameter =
            2 * cell_radius; // diameter of a cell
    const size_t N = 5; // initial number of cells
    double l_filo_y = 27.5;//2; // sensing radius, filopodia + cell radius
    double l_filo_x = 27.5; // sensing radius, it will have to be rescaled when domain grows
    double l_filo_x_in = l_filo_x; // this value is used for rescaling when domain grows based on initial value
    double l_filo_max = 45; // this is the length when two cells which were previously in a chain become dettached
    //double diff_conc = 0.1; // sensing threshold, i.e. how much concentration has to be bigger, so that the cell moves in that direction
    int freq_growth = 1; // determines how frequently domain grows (actually not relevant because it will go every timestep)
    int insertion_freq = 1; // determines how frequently new cells are inserted, regulates the density of population
    double speed_l = 0.14;// 0.05;//1;//0.05; // speed of a leader cell
    double increase_fol_speed = 1.3;
    double speed_f = increase_fol_speed * speed_l;//0.05;//0.1;//0.08; // speed of a follower cell
    double dettach_prob = 0.5; // probability that a follower cell which is on trail looses the trail
    double chemo_leader = 0.9; //0.5; // phenotypic switching happens when the concentration of chemoattractant is higher than this (presentation video 0.95), no phenotypic switching
    double eps = 1; // for phenotypic switching, the distance has to be that much higher
    const int filo_number = 3; // number of filopodia sent
    int same_dir = 0; // number of steps in the same direction +1, because if 0, then only one step in the same direction
    bool random_pers = true; // persistent movement also when the cell moves randomly
    int count_dir = 0; // this is to count the number of times the cell moved the same direction, up to same_dir for each cell
    double lam = 1.0;//72/(100)/10; // to 1000 /h chemoattractant internalisation

    int value = 0; // value of the Gamma(value), where Gamma is close to a cell center

    // distance to the track parameters
    double dist_thres = 1;
    int closest_time;
    int leader_track;


    vdouble2 xposi; // to store positions


    //switch of growth profile
    double first_profile = 26.995; //41.995;//;just once if timestep is 0.01
    double second_profile = 27.0095;//42.0095;//


    // int n_seed = 0;
    // double diff_conc = 0.05;


    // domain growth parameters


    /*
    * strain rate
    * */


    //piecewise constant, two parts
    // 1 part n_faster times faster than the other part
    double n_faster = 2.0;

    double thetasmall = 0.75; // first thetasmall is growing
    int theta1 = int(thetasmall * length_x);

    double alpha1;
    double alpha2;

    //check whether first part grows
    double thetasmalltemp = thetasmall;
    if (first_part_grows ==true){
        double xvar = final_length/ (n_faster * double(length_x)*thetasmalltemp + double(length_x)*(1-thetasmalltemp)); // solve: 2 *xvar * length_x * thetasmall + x * length(1-thetasmall) = final_length

        double ratio1 = n_faster * double(length_x)*thetasmalltemp * xvar / (double(length_x)*thetasmalltemp);

        alpha1 = log (ratio1)/final_time;

        double ratio2 =  double(length_x)*(1- thetasmalltemp) * xvar / (double(length_x)*(1- thetasmalltemp));

        alpha2 = log (ratio2)/final_time;

    }

    else if(first_part_grows == false){
        double xvar = final_length/ (double(length_x)*thetasmall + n_faster*double(length_x)*(1-thetasmall)); // solve: 2 *xvar * length_x * thetasmall + x * length(1-thetasmall) = final_length

        double ratio1 = double(length_x)*thetasmall * xvar / (double(length_x)*thetasmall);

        alpha1 = log (ratio1)/final_time;

        double ratio2 =  double(length_x)*(1- thetasmall) * xvar*n_faster / (double(length_x)*(1- thetasmall));

        alpha2 = log (ratio2)/final_time;
    }




    VectorXd strain = VectorXd::Zero(length_x);


    // constant
//    for (int i = 0; i < length_x; i++) {
//        strain(i) = alpha;// 0;//linear_par * double(theta1) / double(space_grid_controller);//0;
//    }

    //cout << alpha1 << alpha1 << endl;
//    alpha1 = 0.0301;
//    alpha2 = 0.0;

    // first part it is linear growth
    for (int i = 0; i < theta1; i++) {
        strain(i) = alpha1;//linear_par * double(theta1) /
        // double(space_grid_controller);//linear_par * (double(i) / double(space_grid_controller));
    }


    // second part is constant
    for (int i = theta1; i < length_x; i++) {
        strain(i) = alpha2;// 0.002;//0.5 * linear_par * double(theta1) / double(space_grid_controller); // constant to where it was
        //strain(i,j) = linear_par*theta1/(theta1- (theta2-1))*(i-(theta2-1)); // linearly decreasing, if I do this do not forget to change Gamma
    }

    cout << "strain(0) " << strain(0) << endl;
    cout << "strain(1) " << strain(length_x-2) << endl;





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
    VectorXd Gamma_temp = VectorXd::Zero(length_x);     // for switch in growth profiles
    VectorXd Gamma_base = VectorXd::Zero(length_x);     // for switch in growth profiles

    for (int i = 0; i < length_x; i++) {
        Gamma_x(i) = exp(0 * strain(i));
        //Gamma(i) = Gamma_initial * i / space_grid_controller;
        Gamma(i) = i;
        Gamma_old(i) = Gamma(i);
    }

    double Gamma_initial = Gamma(length_x - 1);


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

    for (int i = 0; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            chemo(i, j) = 1;//cos(
            // M_PI * i / space_grid_controller);//1;//C0 - 0.5 * cos( M_PI * i/space_grid_controller * n);
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
//        chemo_3col(i, 0) = Gamma_initial * chemo_3col_ind(i, 0) / double(space_grid_controller);
        chemo_3col(i, 0) = chemo_3col_ind(i, 0);// / double(space_grid_controller);

    }



    // u column
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
    }


    // y coordinates, 1D so nothing changes
    for (int i = 0; i < length_x * length_y; i++) {
        // chemo_3col(i, 1) = y_init * chemo_3col_ind(i, 1) / double(space_grid_controller);
        chemo_3col(i, 1) = chemo_3col_ind(i, 1);// / double(space_grid_controller);
    }


    int counter = 0;

    while (t < final_time) {





        /*
         * After some time, domain growth profile changes!!!
         *
         * */
//
        if (t > first_profile && t < second_profile) {

            bool first_part_grows = false; //false if final part grows faster, change



            n_faster = 2.0;

            thetasmall = 0.5; // first thetasmall is growing
            theta1 = int(thetasmall * length_x);

            alpha1;
            alpha2;

            //check whether first part grows
            double thetasmalltemp = thetasmall;
            if (first_part_grows == true) {
                double xvar = final_length / (n_faster * double(length_x) * thetasmalltemp + double(length_x) * (1 -
                                                                                                                 thetasmalltemp)); // solve: 2 *xvar * length_x * thetasmall + x * length(1-thetasmall) = final_length

                double ratio1 =
                        n_faster * double(length_x) * thetasmalltemp * xvar / (double(length_x) * thetasmalltemp);

                alpha1 = log(ratio1) / final_time;

                double ratio2 =
                        double(length_x) * (1 - thetasmalltemp) * xvar / (double(length_x) * (1 - thetasmalltemp));

                alpha2 = log(ratio2) / final_time;

            }
            else if (first_part_grows == false) {

                double xvar = final_length / (double(length_x) * thetasmall + n_faster * double(length_x) * (1 -
                                                                                                             thetasmall)); // solve: 2 *xvar * length_x * thetasmall + x * length(1-thetasmall) = final_length

                double ratio1 = double(length_x) * thetasmall * xvar / (double(length_x) * thetasmall);

                alpha1 = log(ratio1) / final_time;

                double ratio2 =
                        double(length_x) * (1 - thetasmall) * xvar * n_faster / (double(length_x) * (1 - thetasmall));

                alpha2 = log(ratio2) / final_time;
            }


            // constant
//    for (int i = 0; i < length_x; i++) {
//        strain(i) = alpha;// 0;//linear_par * double(theta1) / double(space_grid_controller);//0;
//    }



            // first part it is linear growth
            for (int i = 0; i < theta1; i++) {
                strain(i) = alpha1;//linear_par * double(theta1) /
                // double(space_grid_controller);//linear_par * (double(i) / double(space_grid_controller));
            }


            // second part is constant
            for (int i = theta1; i < length_x; i++) {
                strain(i) = alpha2;// 0.002;//0.5 * linear_par * double(theta1) / double(space_grid_controller); // constant to where it was
                //strain(i,j) = linear_par*theta1/(theta1- (theta2-1))*(i-(theta2-1)); // linearly decreasing, if I do this do not forget to change Gamma
            }
        }


        /*
 *
 * Domain growth, and update cell positions
 */




        // update the strain rate


        if (t < first_profile) {
            for (int i = 0; i < length_x; i++) {
                Gamma_x(i) = exp(t * strain(i));
            }
        }

//changed here!!!!
        if (t > first_profile) {
            for (int i = 0; i < length_x; i++) {
                Gamma_x(i) = exp((t-first_profile) * strain(i));
            }

        }

        /*
         * arbitrary Gamma function
         * */


        Gamma(0) = 0; // this is assumption, since I cannot calculate it
        //cout << "Gamma(0) " << 0 << " value " << Gamma(0) << endl;


        if (t < first_profile) {
            for (int i = 1; i < length_x; i++) {

                Gamma(i) = Gamma_x(i) * dx + Gamma(i - 1);
            }
            Gamma_base = Gamma;
        }


        Gamma_temp(0) = 0; // this is assumption, since I cannot calculate it

        double total_increase = 0;
        double smallest_increase = 0;
        double total_sum = Gamma_x(0);
        double factor = 0; // sum all gamma over this factor is equal to total change

        if (t > first_profile) {
            // old version which worked
            //this is for general case
//            for (int i = 1; i < length_x;i++){
//
//                Gamma_temp(i) = Gamma_x(i) * dx + Gamma_temp(i-1);
//                total_sum = total_sum +Gamma_x(i);
//            }
//            total_increase = Gamma_temp(length_x-1) - Gamma(length_x-1);
//
//            factor = total_sum/total_increase;
//
//
//            for (int i = 1; i < length_x;i++){
//                Gamma(i) = Gamma_old(i) + Gamma_temp(i)/factor;
//            }

// new version
            for (int i = 1; i < theta1;i++){
                Gamma(i) = Gamma_base(i)* Gamma_x(i);
            }
            for (int i = theta1; i < length_x;i++){

                Gamma(i) =  Gamma_base(i)* Gamma_x(i) + Gamma_base(theta1)*(Gamma_x(1)-Gamma_x(i));
            }


        }

        // update time step;

        t = t + dt;

        counter = counter + 1;




        /*
         * Domain growth
         * */

//        Lt = thetasmall + thetasmall * exp(alpha1 * t);
//        Lt = Gamma(length_x - 1);
//
//
//        //   Ltdot = alpha * thetasmall * exp(alpha * t); // piecewise, exact
//
//        Ltdot = (Lt - Lt_old) / dt; //could be used for both, especially should be used if derivative is unknown
//        //    cout << "diff dot " << Ltdot << endl;
//
//        Lt_old = Lt;


        // I need Gamma_t for cos verification as well

        for (int i = 0; i < length_x; ++i) {

            Gamma_t(i) = (Gamma(i) - Gamma_old(i)) / dt;

        }



        /// update positions uniformly based on the domain growth

        vdouble2 x; // use variable x for the position of cells
        double x0 = 0;
        int pos;




        Gamma_old = Gamma;




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






        // u column
        for (int i = 0; i < length_x * length_y; i++) {
            chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
        }

        if (counter % 100 == 0) {






            //ofstream output("matrix_FIRST_025theta" + to_string(int(round(t))) + ".csv");
            ofstream output("Domain075first05finalV2" + to_string(t) + ".csv");


            output << "x, y, z, u" << "\n" << endl;



            //output << "x, y, z, u" << "\n" << endl;


            for (int i = 0; i < length_x * length_y; i++) {
                for (int j = 0; j < 4; j++) {
                    output << chemo_3col(i, j) << ", ";
                }
                output << "\n" << endl;
            }



        }



    }




    }




/*
 * main for proportions in different sections
 */


// parameter analysis
int main() {

    const int number_parameters = 1; // parameter range
    const int sim_num = 1;

    //VectorXd store_chains;
    VectorXi vector_check_length = proportions(0.05, 0); //just to know what the length is

    //int num_parts = vector_check_length.size(); // number of parts that I partition my domain
    //cout << "length " << vector_check_length.size() << endl;
    int num_parts = 18; // since 18 is domain partition for time 54 and 1 for chains
    MatrixXf sum_of_all = MatrixXf::Zero(num_parts, number_parameters); // sum of the values over all simulations

    //initialise the matrix to store the values
    MatrixXi numbers = MatrixXi::Zero(num_parts, number_parameters);

    // n would correspond to different seeds
    // parallel programming




}
