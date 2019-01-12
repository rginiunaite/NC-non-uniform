//
// Created by rasa on 07/01/19.
//

/*
 * Reaction-diffusion equation on a non-uniformly growing domain, 2D, explicit method, piece-wise non-uniform growth, piecewise two parts, piecewise three parts with one linear
 * Cell dynamics included
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




//VectorXi proportions(double diff_conc, int n_seed) {
int  main(){


// parameters

// model parameters
// model parameters

/*
 * Initial length of the domain is 300 \mu m, I will create a matrix of chemoattractant which will be of length 300/10= 30.
 * All the length parameters will be scaled consistently, i.e. /10.
 */

    double space_grid_controller = 100.0;

    double domain_length = 3.0; //this variable is for the actual domain length, since it will be increasing
    double Lt_old = domain_length;
    int length_x =
            int(domain_length) * int(space_grid_controller); // length in x direction of the chemoattractant matrix
    double initial_domain_length = domain_length;
    const int length_y = int(1.2*space_grid_controller); // length in y direction of the chemoattractant matrix
    const double final_time = 72; // number of timesteps, 1min - 1timestep, from 6h tp 24hours.

// parameters for the dynamics of chemoattractant concentration

    double D = 0.00001;//0.00001;//0.05; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0.0; // initialise time
    double dt = 0.05;//before was 0.05 // time step
    double dt_init = dt;
    int number_time = int(1 / dt_init); // how many timesteps in 1min, which is the actual simulation timestep
    double dx = 1.0/space_grid_controller; // space step in x direction, double to be consistent with other types
    double dy = 1.0/space_grid_controller; // space step in y direction
    double kai = 0.00001;//0;//0.1 // to 1 /h production rate of chemoattractant
//    double y_init = 120.0;

    // reaction rate
    double k_reac = 0;//0.00001;//0.1;//0.105;//0.03;//.205; // reaction term


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
    double speed_l = 0.5;// 0.05;//1;//0.05; // speed of a leader cell
    double increase_fol_speed = 1.3;
    double speed_f = increase_fol_speed * speed_l;//0.05;//0.1;//0.08; // speed of a follower cell
    double dettach_prob = 0.5; // probability that a follower cell which is on trail looses the trail
    double chemo_leader = 0.9; //0.5; // phenotypic switching happens when the concentration of chemoattractant is higher than this (presentation video 0.95), no phenotypic switching
    double eps = 1; // for phenotypic switching, the distance has to be that much higher
    const int filo_number = 3; // number of filopodia sent
    int same_dir = 0; // number of steps in the same direction +1, because if 0, then only one step in the same direction
    bool random_pers = true; // persistent movement also when the cell moves randomly
    int count_dir = 0; // this is to count the number of times the cell moved the same direction, up to same_dir for each cell
    double lam = 2.0;//72/(100)/10; // to 1000 /h chemoattractant internalisation


    // distance to the track parameters
    double dist_thres = 1;
    int closest_time;
    int leader_track;



    int n_seed = 0;
    double diff_conc = 0.05;


    // domain growth parameters


    /*
    * strain rate
    * */

    //double linear_par = 0.0001;//05;


    // for comparison with analytical
    double alpha = 0.02563;//for 0.05 time step and 72 final time, theta1 = 0.5
    //double alpha = 0.03412; // for 0.05 time step, theta1 = 0.75
    //double alpha = 1.706;// dt = 0.001, 0.75


    VectorXd strain = VectorXd::Zero(length_x);


    // constant
//    for (int i = 0; i < length_x; i++) {
//        strain(i) = alpha;// 0;//linear_par * double(theta1) / double(space_grid_controller);//0;
//    }

    //piecewise constant, two parts

    int theta1 = int(0.5 * length_x);
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
        //Gamma(i) = Gamma_initial * i / space_grid_controller;
        Gamma(i) = i;
        Gamma_old(i) = Gamma(i);
    }

    double Gamma_initial = Gamma(length_x-1);


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


    // save data to plot chemoattractant concentration
    ofstream outputinit("matrix_non_uniform0.csv");

    outputinit << "x, y, z, u" << "\n" << endl;


    for (int i = 0; i < length_x * length_y; i++) {
        for (int j = 0; j < 4; j++) {
            outputinit << chemo_3col(i, j) << ", ";
        }
        outputinit << "\n" << endl;
    }


    //ofstream output2("track_point" + to_string(t) + ".csv");

//    ofstream output2in("track_point0.csv");
//
//    output2in << Gamma(length_x / 2) << endl;
//
//    ofstream output3in("track2_point0.csv");
//
//    output3in << Gamma(int(length_x / 4)) << endl;
//
//    ofstream output4in("track3_point.csv");
//
//    output4in << Gamma(3 * int(length_x / 4)) << endl;
//


    /*
     * 2D domain with a few randomly placed particles
     */

    /*
     * initial cells of fixed radius
     */

    //ABORIA_VARIABLE(velocity,vdouble2,"velocity")
    ABORIA_VARIABLE(radius, double, "radius")
    ABORIA_VARIABLE(direction, vdouble2, "direction")// stores the direction a particle moved
    ABORIA_VARIABLE(persistence_extent, int,
                    "persistence extent")// stores whether cell moves only one step in current direction or in a process of moving persistently
    ABORIA_VARIABLE(same_dir_step, int,
                    "same dir step")// the number which stores how many steps in the same direction are made.
    ABORIA_VARIABLE(attached_to_id, int, "attached_to_id")
    ABORIA_VARIABLE(type, int, "type") // 0 if a cell is a leader, 1 if follower
    ABORIA_VARIABLE(chain_type, int, "chain_type") // leaders form different chain types
    ABORIA_VARIABLE(chain, int, "chain") // stores whether a follower is part of the chain or no, 0 if it is not part of
    // the chain and then increasing integer if it is. If it is attached to a leader, it is 1, and then increasing order.
    // stores the distance to the closest neighbour, if less than thresold
    typedef Particles<std::tuple<radius, type, attached_to_id, direction, chain, chain_type, persistence_extent, same_dir_step>, 2> particle_type; // 2 stands for dimension

    // will use stored value of the position of a particle
    typedef particle_type::position position;

    // initialise the number of particles
    particle_type particles(2*N);

    // initialise random number generator for particles entering the domain, appearing at the start in x and uniformly in y
    std::default_random_engine gen;
    std::uniform_real_distribution<double> uniform(cell_radius, length_y - cell_radius);


    /*
     * compact initialisation of particles
     */

    for (int i = 0; i < N; ++i) {


        get<radius>(particles[i]) = cell_radius;
        get<type>(particles[i]) = 0; // initially all cells are leaders

        //get<position>(p) = vdouble2(cell_radius,(i+1)*diameter); // x=2, uniformly in y
        get<position>(particles[i]) = vdouble2(cell_radius + theta1, (i + 1) * double(length_y - 1) / double(N) -
                                                            0.5 * double(length_y - 1) /
                                                            double(N)); // x=2, uniformly in y
        get<persistence_extent>(particles[i]) = 0;
        get<same_dir_step>(particles[i]) = 0;

    }

    for (int i = 0; i < N; ++i) {


        get<radius>(particles[i+N]) = cell_radius;
        get<type>(particles[i+N]) = 0; // initially all cells are leaders

        //get<position>(p) = vdouble2(cell_radius,(i+1)*diameter); // x=2, uniformly in y
        get<position>(particles[i+N]) = vdouble2(theta1-cell_radius, (i + 1) * double(length_y - 1) / double(N) -
                                                                     0.5 * double(length_y - 1) /
                                                                     double(N)); // x=2, uniformly in y
        get<persistence_extent>(particles[i+N]) = 0;
        get<same_dir_step>(particles[i+N]) = 0;

    }

    // initialise neighbourhood search, note that the domain will grow in x direction, so I initialise larger domain
    particles.init_neighbour_search(vdouble2(0, 0), 5 * vdouble2(length_x, length_y), vbool2(false, false));

    // save particles before they move

    vtkWriteGrid("particles", t, particles.get_grid(true));

    // initialise random number generator to obtain random number between 0 and 2*pi
    std::default_random_engine gen1;
    gen1.seed(t * n_seed); // choose different seeds to obtain different random numbers
    std::uniform_real_distribution<double> uniformpi(0, 2 * M_PI);












    //for each timestep
    while (t < final_time) {



//              insert new cells


//            bool free_position = false;
//            particle_type::value_type f;
//            //get<radius>(f) = cell_radius;
//
//
//            get<position>(f) = vdouble2(cell_radius, uniform(gen)); // x=2, uniformly in y
//            free_position = true;
//            /*
//             * loop over all neighbouring leaders within "dem_diameter" distance
//             */
//            for (auto tpl = euclidean_search(particles.get_query(), get<position>(f), diameter); tpl != false; ++tpl) {
//
//                vdouble2 diffx = tpl.dx();
//
//                if (diffx.norm() < diameter) {
//                    free_position = false;
//                    break;
//                }
//            }
//
//            // our assumption that all new cells are followers
//            get<type>(f) = 0;
//
//
//            if (free_position) {
//                get<chain>(f) = 0;
//                get<chain_type>(f) = -1;
//                get<attached_to_id>(f) = -1;
//                particles.push_back(f);
//            }
//
//
//        particles.update_positions();



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
//            Gamma(i) =
//                    Gamma_initial * (double(i) / double(space_grid_controller)) * (Gamma_x(i)); // linearly increasing
            Gamma(i) = (double(i)) * (Gamma_x(i)); // linearly increasing

        }



        // second part is constant
        for (int i = theta1; i < length_x; i++) {
//
//            Gamma(i) = Gamma(theta1 - 1) +
//                       Gamma_initial * (double(i) / double(space_grid_controller) - double(theta1 - 1) /
//                                                                                    double(space_grid_controller)) *
//                       Gamma_x(i);

            Gamma(i) = Gamma(theta1 - 1) + (double(i)- double(theta1 - 1)) * Gamma_x(i);
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



        /// update positions uniformly based on the domain growth

        vdouble2 x; // use variable x for the position of cells
        int pos;

        for (int i = 0; i < particles.size(); i++) {

            x = get<position>(particles[i]);
// wrong!!!!
//            for(Index j=0; j<Gamma_old.size(); ++j){
//                if (x[0] > Gamma_old(j)){
//                    pos = j;
//                }
//                else{
//                    break;
//                }
//            }
//
//            cout << "x pos " << x[0] << endl;
//            cout << " gamma_old(j) " << Gamma_old(pos) << endl;
//            cout << "gamma(j) " << Gamma(pos) << endl;

// since I do not know how to do it for general case, I will do it for my specific

        if (x[0] > theta1){
            get<position>(particles)[i] *= vdouble2((Gamma(length_x-1)- theta1)/(Gamma_old(length_x-1) - theta1) , 1); // update position based on changes in Gamma
        }

        }


        Gamma_old = Gamma;




        // internalisation rate



        // internalisation
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                //go through all the cells
                for (int k = 0; k < particles.size(); k++) {
                    // leaders
                    //for (int k = 0; k < N; k++) {
                    vdouble2 x;
                    x = get<position>(particles[k]);
                    intern(i, j) = intern(i, j) + exp(-((Gamma(i) - x[0]) *
                                                        (Gamma(i) - x[0]) +
                                                        (j - x[1]) * (j - x[1])) /
                                                      (2 * cell_radius * cell_radius)); // mapping to fixed domain
                }
            }
        }



        // inner coefficients


        for (int i = 1; i < length_x - 1; ++i) {
            for (int j = 1; j < length_y - 1; ++j) {

                chemo_new(i, j) = dt * (D * 1.0 / (2.0 * dx * dx * Gamma_x(i)) *
                                        ((1.0 / Gamma_x(i) + 1.0 / Gamma_x(i + 1)) * (chemo(i + 1) - chemo(i)) -
                                         (chemo(i, j) - chemo(i - 1, j)) * (1.0 / Gamma_x(i) + 1.0 / Gamma_x(i - 1))) +
                                        D * (chemo(i, j + 1) - 2 * chemo(i, j) + chemo(i, j - 1)) / (dy * dy) -
                                        (chemo(i, j) * lam / (2 * M_PI * cell_radius * cell_radius)) * intern(i, j) +
                                        chemo(i, j) * k_reac - strain(i) * chemo(i, j)) + chemo(i, j);
            }
        }
        // i = 0, will have to loop over all y when I will move to two dimensions



        for (int i = 0; i < length_y; i++) {
            chemo_new(0, i) = chemo_new(1, i);
            chemo_new(length_x - 1, i) = chemo_new(length_x - 2, i);

        }

        for (int i = 0; i < length_x; i++) {
            chemo_new(i, 0) = chemo_new(i, 1);
            chemo_new(i, length_y - 1) = chemo_new(i, length_y - 2);
        }


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





        /*
         * Update the position of particles
         * */


        //  create a random list of cell ids
        int check_rep = 0; // check for repetitions, 0 no rep, 1 rep


        std::default_random_engine gen2;
        gen2.seed(t * n_seed); // different seeds
        std::uniform_real_distribution<double> uniform_particles(0, particles.size()); // can only move forward

        VectorXi particle_id = VectorXi::Zero(particles.size());

        for (int i = 0; i < particles.size(); i++) {

            check_rep = 1; // set to 1 to enter the while loop
            while (check_rep == 1) {
                check_rep = 0; // it is initially zero and then will be changed to 1 if it is equivalent to others
                particle_id(i) = uniform_particles(gen2);


                for (int j = 0; j < i; j++) {
                    if (particle_id(i) == particle_id(j)) { check_rep = 1; }
                }
            }
        }

        // update the position of all particles in a random order created above

        for (int j = 0; j < particles.size(); j++) {

            /*
            * phenotypic switching, based on chemoattractant concentration in front, +0.5
            */
//
//            vdouble2 coord = get<position>(particles[particle_id(j)]);
//
//            // rescaled coord
//
//            double rescaled_coord;
//
//            rescaled_coord = (length_x / domain_length)*coord[0];
//
//            double chemo_in_front = chemo(round(rescaled_coord), round(coord[1]));
//            //cout << "chemo in front " << old_chemo << endl;
//
//
//            // if high concentration cells become leaders
//            if (chemo_in_front > chemo_leader ){
//                get<type>(particles[particle_id(j)]) = 0;
//            }
//            else{
//                get<type>(particles[particle_id(j)]) = 1;
//            }


            // if a particle is a leader


            vdouble2 x; // use variable x for the position of cells
            x = get<position>(particles[particle_id(j)]);

            if (get<type>(particles[particle_id(j)]) == 0) {

                vdouble2 x; // use variable x for the position of cells
                x = get<position>(particles[particle_id(j)]);


                double x_in; // x coordinate in initial domain length scale
                if (x[0] < theta1){
                    x_in = x[0];
                }
                else{
                    x_in = x[0]* (Gamma_initial -theta1)/(Gamma(length_x-1)- theta1) ;
                    l_filo_x = l_filo_x_in * (Gamma_initial -theta1)/(Gamma(length_x-1)- theta1) ; // rescale the length of filopodia as well
                }




                // if it is still in the process of moving in the same direction
                if (get<persistence_extent>(particles[particle_id(j)]) == 1) {


                    x += get<direction>(particles[particle_id(j)]);

                    if (x[0] < theta1){
                        x_in = x[0];
                    }
                    else{
                        x_in = x[0]* (Gamma_initial -theta1)/(Gamma(length_x-1)- theta1) ;
                    }


                    bool free_position = true; // check if the neighbouring position is free

                    // check if there are other particles in the position where the particle wants to move
                    for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {
                        if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            free_position = false;
                        }
                    }

                    // check that the position they want to move to is free and not out of bounds
                    if (free_position && x[0] > cell_radius && x[0] < Gamma(length_x - 1) && (x[1]) > cell_radius &&
                        (x[1]) < length_y - cell_radius) {
                        // if that is the case, move into that position
                        get<position>(particles)[particle_id(j)] +=
                                get<direction>(particles)[particle_id(j)];
                    }
                    get<same_dir_step>(particles)[particle_id(
                            j)] += 1; // add regardless whether the step happened or no to that count of the number of movement in the same direction

                }


                // if a particle is not in a sequence of persistent steps
                if (get<persistence_extent>(particles[particle_id(j)]) == 0) {



                    // create an array to store random directions
                    array<double, filo_number + 1> random_angle;

                    // choose the number of angles where the filopodia is sent
                    for (int k = 0; k < filo_number + 1; k++) {

                        double random_angle_tem = uniformpi(gen1);
                        int sign_x_tem, sign_y_tem;

                        random_angle_tem = uniformpi(gen1);
                        random_angle[k] = random_angle_tem;

                    }


                    // choose which direction to move


                    // store variables for concentration at new locations


                    double old_chemo = chemo((round(x_in)), round(x)[1]);
                    array<double, filo_number> new_chemo;


                    for (int i = 0; i < filo_number; i++) {

                        if (round((x_in + sin(random_angle[i]) * l_filo_x)) < 0 ||
                            round((x_in + sin(random_angle[i]) * l_filo_x)) >
                            length_x - 1 || round(x[1] + cos(random_angle[i]) * l_filo_y) < 0 ||
                            round(x[1] + cos(random_angle[i]) * l_filo_y) > length_y - 1) {
                            new_chemo[i] = 0;
                        } else {
                            new_chemo[i] = chemo(round((x_in + sin(random_angle[i]) * l_filo_x)),
                                                 round(x[1] + cos(random_angle[i]) * l_filo_y));
                        }

                    }


                    // find maximum concentration of chemoattractant

                    int chemo_max_number = 0;

                    for (int i = 1; i < filo_number; i++) {
                        if (new_chemo[chemo_max_number] < new_chemo[i]) {
                            chemo_max_number = i;
                        }
                    }

                    // if the concentration in a new place is relatively higher than the old one (diff_conc determines that threshold), move that way
                    if ((new_chemo[chemo_max_number] - old_chemo) / sqrt(old_chemo) > diff_conc) {

                        count_dir += 1;

                        x += speed_l *
                             vdouble2(sin(random_angle[chemo_max_number]), cos(random_angle[chemo_max_number]));



                        if (x[0] < theta1){
                            x_in = x[0];
                        }
                        else{
                            x_in = x[0]* (Gamma_initial -theta1)/(Gamma(length_x-1)- theta1) ;
                        }

                        bool free_position = true; // check if the neighbouring position is free

                        // check if the position the particle wants to move is free
                        for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                            if (get<id>(*k) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                free_position = false;
                            }
                        }


                        // if the position they want to move to is free and not out of bounds, move that direction
                        if (free_position && x[0] > cell_radius && x[0] < Gamma(length_x - 1) && (x[1]) > cell_radius &&
                            (x[1]) < length_y - cell_radius) {
                            get<position>(particles)[particle_id(j)] +=
                                    speed_l * vdouble2(sin(random_angle[chemo_max_number]),
                                                       cos(random_angle[chemo_max_number])); // update if nothing is in the next position
                            get<direction>(particles)[particle_id(j)] =
                                    speed_l * vdouble2(sin(random_angle[chemo_max_number]),
                                                       cos(random_angle[chemo_max_number]));

                            // if there is some kind of tendency to move persistently
                            if (same_dir > 0) {
                                get<persistence_extent>(particles[particle_id(
                                        j)]) = 1; // assume for now that it also becomes peristent in random direction

                            }

                        }


                    }

                        // if the concentration is not higher, move in random direction
                    else {


                        x += speed_l * vdouble2(sin(random_angle[filo_number]), cos(random_angle[filo_number]));



                        if (x[0] < theta1){
                            x_in = x[0];
                        }
                        else{
                            x_in = x[0]* (Gamma_initial -theta1)/(Gamma(length_x-1)- theta1) ;
                        }


                        bool free_position = true; // check if the neighbouring position is free

                        // if this loop is entered, it means that there is another cell where I want to move
                        for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                            if (get<id>(*k) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                free_position = false;
                            }
                        }


                        // update the position if the place they want to move to is free and not out of bounds
                        if (free_position && x[0] > cell_radius && x[0] < Gamma(length_x - 1) && (x[1]) > cell_radius &&
                            (x[1]) < length_y - cell_radius) {
                            get<position>(particles)[particle_id(j)] +=
                                    speed_l * vdouble2(sin(random_angle[filo_number]),
                                                       cos(random_angle[filo_number])); // update if nothing is in the next position
                            get<direction>(particles)[particle_id(j)] =
                                    speed_l * vdouble2(sin(random_angle[filo_number]),
                                                       cos(random_angle[filo_number]));
                            // if particles start moving persistently in all directions
                            if (random_pers) {
                                if (same_dir > 0) {
                                    get<persistence_extent>(particles[particle_id(
                                            j)]) = 1; // assume for now that it also becomes peristent in random direction

                                }
                            }

                        }

                    }

                }

                // check if it is not the end of moving in the same direction
                if (get<same_dir_step>(particles)[particle_id(j)] > same_dir) {
                    get<persistence_extent>(particles)[particle_id(j)] = 0;
                    get<same_dir_step>(particles[particle_id(j)]) = 0;
                }

            }




            // if a particle is a follower
            if (get<type>(particles[particle_id(j)]) == 1) {

                vdouble2 x;
                x = get<position>(particles[particle_id(j)]);

                double x_in; // x coordinate in initial domain length scale


                if (x[0] < theta1){
                    x_in = x[0];
                }
                else{
                    x_in = x[0]* (Gamma_initial -theta1)/(Gamma(length_x-1)- theta1) ;
                }
                // if the particle is part of the chain
                if (get<chain>(particles[particle_id(j)]) > 0) {


                    // check if it is not too far from the cell it was following

                    vdouble2 dist;

                    dist = get<position>(particles[particle_id(j)]) -
                           get<position>(particles[get<attached_to_id>(particles[particle_id(j)])]);

                    // if it is sufficiently far dettach the cell
                    if (dist.norm() > l_filo_max) {
                        get<chain>(particles[particle_id(j)]) = 0;
                        //dettach also all the cells that are behind it, so that other cells would not be attached to this chain
                        for (int i = 0; i < particles.size(); ++i) {
                            if (get<chain_type>(particles[i]) == get<chain_type>(particles)[particle_id(j)]) {
                                // get<chain_type>(particles)[i] = -1;
                                get<chain>(particles[i]) = 0;
                            }

                        }
                        // get<chain_type>(particles)[particle_id(j)] = -1;
                    }


                    // direction the same as of the cell it is attached to
                    get<direction>(particles)[particle_id(j)] = get<direction>(particles)[get<attached_to_id>(
                            particles[particle_id(j)])];

                    //try to move in the same direction as the cell it is attached to
                    vdouble2 x_chain = x + increase_fol_speed * get<direction>(particles)[particle_id(j)];

                    double x_in_chain; // scaled coordinate

                    if (x_chain[0] < theta1){
                        x_in_chain = x_chain[0];
                    }
                    else{
                        x_in_chain = x_chain[0]* (Gamma_initial -theta1)/(Gamma(length_x-1)- theta1) ;
                    }


                    bool free_position = true;

                    // check if the position it wants to move to is free
                    for (auto pos = euclidean_search(particles.get_query(), x_chain, diameter);
                         pos != false; ++pos) {
                        if (get<id>(*pos) !=
                            get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            free_position = false;
                        }
                    }


                    // update the position if the place they want to move to is free and not out of bounds
                    if (free_position && x_chain[0] > cell_radius && x_chain[0] < Gamma(length_x - 1) && (x_chain[1]) > cell_radius &&
                        (x_chain[1]) < length_y - cell_radius) {
                        get<position>(particles)[particle_id(j)] +=
                                increase_fol_speed * get<direction>(particles[particle_id(j)]);

                    }

//                        /*
//                           * NOT IN THE SUMMARY I SENT THEM
//                           */
//                    if (free_position == false){
//                        get<chain>(particles[particle_id(j)]) = 0;
//                    }





                }

                // if the cell is not part of the chain
                if (get<chain>(particles[particle_id(j)]) == 0) {


                    /* check if there are any cells distance l_filo_y apart
                    * it can be either a leader or a follower already in a chain
                    */


                    // try to find a close leader
                    for (auto k = euclidean_search(particles.get_query(), x, l_filo_x_in); k != false; ++k) {


                        if (get<type>(*k) == 0) { // if it is close to a leader
                            get<direction>(particles)[particle_id(j)] = get<direction>(*k); // set the same direction
                            get<chain>(particles)[particle_id(j)] = 1; // note that it is directly attached to a leader
                            get<attached_to_id>(particles)[particle_id(j)] = get<id>(
                                    *k); // note the id of the particle it is attached to
                            get<chain_type>(particles)[particle_id(j)] = get<id>(
                                    *k); // chain type is the id of the leader
                        }

                    }


                    // if it hasn't found a leader nearby,
                    // try to find a close follower which is in a chain contact with a leader

                    if (get<chain>(particles)[particle_id(j)] != 1) {
                        for (auto k = euclidean_search(particles.get_query(), x, l_filo_y); k != false; ++k) {

                            // if it is close to a follower that is part of the chain
                            if (get<type>(*k) == 1 && get<chain>(*k) > 0) {

                                if (get<id>(*k) != get<id>(particles[particle_id(j)])) {
                                    //check if there is a leader in front of the chain
                                    //for (int i= 0; i< particles.size(); i++){//go through all the particles
                                    //  if (get<chain_type>(particles)[i] == get<chain_type>(*k) && get<chain>(particles)[i] == 1){// if they are in the same chain and directly attached to a leader

                                    get<direction>(particles)[particle_id(j)] = get<direction>(*k);
                                    get<chain>(particles)[particle_id(j)] =
                                            get<chain>(*k) + 1; // it is subsequent member of the chain
                                    get<attached_to_id>(particles)[particle_id(j)] = get<id>(
                                            *k); // id of the particle it is attached to
                                    get<chain_type>(particles)[particle_id(j)] = get<chain_type>(*k); // chain type is
                                    // the same as the one of the particle it is attached to
                                    //}

                                    //}

                                }


                            }

                        }
                    }

                    // try to move if it has found something

                    if (get<chain>(particles[particle_id(j)]) > 0) {

                        //try to move in the same direction as the cell it is attached to
                        vdouble2 x_chain = x + increase_fol_speed * get<direction>(particles)[particle_id(j)];

                        // Non-uniform domain growth
                        double x_in_chain;

                        if (x_chain[0] < theta1){
                            x_in_chain = x_chain[0];
                        }
                        else{
                            x_in_chain = x_chain[0]* (Gamma_initial -theta1)/(Gamma(length_x-1)- theta1) ;
                        }

                        bool free_position = true;


                        // check if the position it wants to move is free
                        for (auto pos = euclidean_search(particles.get_query(), x_chain, diameter);
                             pos != false; ++pos) {

                            if (get<id>(*pos) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                free_position = false;
                            }
                        }


                        // if the position is free and not out of bounds, move that direction
                        if (free_position &&
                            x_chain[0] > cell_radius &&
                           x_chain[0] < length_x - 1 && (x_chain[1]) > cell_radius &&
                            (x_chain[1]) < length_y - cell_radius) {
                            //cout << "direction " << get<direction>(particles[particle_id(j)]) << endl;
                            get<position>(particles)[particle_id(j)] +=
                                    increase_fol_speed * get<direction>(particles[particle_id(j)]);

                        }
                    }


                    /*
                        * NOT IN THE SUMMARY I SENT THEM
                     */
//                    if (free_position == false){
//                        get<chain>(particles[particle_id(j)]) = 0;
//                    }



                    // if it hasn't found anything close, move randomly

                    if (get<chain>(particles[particle_id(j)]) == 0) {

                        double random_angle = uniformpi(gen1);

                        while (((x_in + sin(random_angle) * l_filo_y) < 0 ||
                                ((x_in + sin(random_angle) * l_filo_y)) >
                                length_x - 1 || (x[1] + cos(random_angle) * l_filo_y) < 0 ||
                                (x[1] + cos(random_angle) * l_filo_y) > length_y - 1)) {
                            random_angle = uniformpi(gen1);
                        }

                        x += speed_f * vdouble2(sin(random_angle), cos(random_angle));


                        if (x[0] < theta1){
                            x_in = x[0];
                        }
                        else{
                            x_in = x[0]* (Gamma_initial -theta1)/(Gamma(length_x-1)- theta1) ;
                        }


                        bool free_position = true; // check if the neighbouring position is free

                        // check if the position the cells want to move to is free
                        for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                            if (get<id>(*k) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                free_position = false;
                            }
                        }

                        // if the position they want to move to is free and not out of bounds, move to that position
                        if (free_position && x[0] > cell_radius && x[0] < Gamma(length_x - 1) && (x[1]) > cell_radius &&
                            (x[1]) < length_y - cell_radius) {
                            get<position>(particles)[particle_id(j)] += speed_f * vdouble2(sin(random_angle),
                                                                                           cos(random_angle)); // update if nothing is in the next position
                            get<direction>(particles)[particle_id(j)] = speed_f * vdouble2(sin(random_angle),
                                                                                           cos(random_angle)); // update direction as well
                        }

                    }

                }

                /* CHECK IF A FOLLOWER DOES NOT BECOME A LEADER
                * Alternative phenotypic switching if a follower overtakes a leader it becomes a leader and that leader follower.
                * I will have to be careful when there will be channels because I will have to choose the closest leader
                * */


                // find the closest leader


                // so that I would not go through all the cells I will choose the ones that are closer to the front

                // minimum position in x of the leaders

                int min_index = 0;

                for (int i = 1; i < N; ++i) {
                    if (get<position>(particles[i])[0] < get<position>(particles[min_index])[0]) {
                        min_index = i;
                    }

                }

                // if a follower is eps further in front than the leader, swap their types
                if (get<position>(particles[particle_id(j)])[0] > get<position>(particles[min_index])[0] + eps) {
                    // find distance to all the leaders
                    double distances[N];
                    vdouble2 dist_vector;
                    //check which one is the closest
                    for (int i = 0; i < N; ++i) {
                        dist_vector = get<position>(particles[particle_id(j)]) - get<position>(particles[i]);
                        distances[i] = dist_vector.norm();

                        int winning_index = 0;
                        for (int i = 1; i < N; ++i) {
                            if (distances[i] < distances[winning_index]) {
                                winning_index = i;
                            }
                        }

                        // if this closest leader is behind that follower, swap them
                        if (get<position>(particles[particle_id(j)])[0] >
                            get<position>(particles[winning_index])[0] + eps) {
                            particle_type::value_type tmp = particles[winning_index];


                            // their position swap

                            vdouble2 temp = get<position>(particles[winning_index]);
                            get<position>(particles[winning_index]) = get<position>(particles[particle_id(j)]);
                            get<position>(particles[particle_id(j)]) = temp;


                        }

                    }


                }

            }

        }

        // update positions
        particles.update_positions();





















        // update chemoattractant profile after the domain grew
//        cout << "t " << t << endl;
//        double ine = t;
//        if(ine == t){
//            cout << "is it in here "<< ine << endl;
//
//        }

        if (counter % 10 == 0) {

            //if (t == 1 || t == 10 || t == 20 ) {
            //cout << "heeere " << endl;
            // save data to plot chemoattractant concentration
            //ofstream output("matrix_non_uniform" + to_string(t) + ".csv");

            // output << "x, y, z, u" << "\n" << endl;

            // save data to plot chemoattractant concentration

            // save at every time step
            #ifdef HAVE_VTK
                        vtkWriteGrid("particlesQuat", t, particles.get_grid(true));
            #endif



            ofstream output("matrix_quater_non_uniform" + to_string(int(t)) + ".csv");


            output << "x, y, z, u" << "\n" << endl;



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

//    /*
// * return the density of cells in domain_partition parts of the domain
// */
//    const int domain_partition = int(domain_length / double(5));; // number of intervalas of 50 \mu m
//
//    VectorXi proportions = VectorXi::Zero(domain_partition); // integer with number of cells in particular part
//
//    double one_part = domain_length / double(domain_partition);
//
//
//    for (int i = 0; i < domain_partition; i++) {
//
//        for (int j = 0; j < particles.size(); j++) {
//            vdouble2 x = get<position>(particles[j]);
//            if (i * one_part < x[0] && x[0] < (i + 1) * one_part) {
//                proportions(i) += 1;
//            }
//        }
//
//    }
//
//    return proportions;


}




/*
 * main for proportions in different sections
 */


//// parameter analysis
//int main() {
//
//    const int number_parameters = 1; // parameter range
//    const int sim_num = 1;
//
//    //VectorXi vector_check_length = proportions(0.05, 2); //just to know what the length is
//    //cout << "prop " << vector_check_length << endl;
//    //int num_parts = vector_check_length.size(); // number of parts that I partition my domain
//    int num_parts = 21; // for 1800 timesteps
//    MatrixXf sum_of_all = MatrixXf::Zero(num_parts, number_parameters); // sum of the values over all simulations
//
//    // n would correspond to different seeds
//    // parallel programming
//#pragma omp parallel for
//    for (int n = 0; n < sim_num; n++) {
//
//        // define parameters that I will change
//
//        array<double, number_parameters> threshold;
//        array<double, 1> slope;
//
//        // set the parameters
//        for (int i = 0; i < number_parameters; i++) {
//            threshold[i] = 0.1;
//            //threshold[i] = 0.005 * (i + 1);// 0.01;
//            //cout << "slope " << slope[i] << endl;
//
//        }
//
//        //initialise the matrix to store the values
//        MatrixXi numbers = MatrixXi::Zero(num_parts, number_parameters);
//
//        //#pragma omp parallel for
//        //        for (int i = 0; i < number_parameters; i++) {
//
//        //for (int j = 0; j < 1; j++) {
//
//        numbers.block(0, 0, num_parts, 1) = proportions(threshold[0], 0);
//
//        //}
//        // }
//
//
//        // This is what I am using for MATLAB
//        ofstream output2("control_case.csv");
//
//        for (int i = 0; i < numbers.rows(); i++) {
//
//            for (int j = 0; j < numbers.cols(); j++) {
//
//                output2 << numbers(i, j) << ", ";
//
//                sum_of_all(i, j) += numbers(i, j);
//
//            }
//            output2 << "\n" << endl;
//        }
//
//    }
//    /*
//    * will store everything in one matrix, the entries will be summed over all simulations
//    */
//
//    ofstream output3("aqp_M9.csv");
//
//    for (int i = 0; i < num_parts; i++) {
//
//        for (int j = 0; j < number_parameters; j++) {
//
//            output3 << sum_of_all(i, j) << ", ";
//
//        }
//        output3 << "\n" << endl;
//    }
//
//
//}
