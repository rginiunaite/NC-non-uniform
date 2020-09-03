//
// Created by giniunaite on 22/07/20.
// non uniform growth, two parts, cell-induced domain growth
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




VectorXi proportions(double diff_conc, int n_seed, double domain_growth_par) {
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
    //const double final_time = 24; //chemical ablation of growth 12hours
    double final_length = 1014;//real 1014
    //double final_length = 520;// chemical ablation of growth, 520 control, 16hours // 400 uniform ablation, 510 linear ablation. For D3 I will try 550 for linear growth


// parameters for the dynamics of chemoattractant concentration

    double D = 2.0;//0.00001;//0.05; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0.0; // initialise time
    double dt = 0.01; // time step
    double dt_init = dt;
    int number_time = int(1.0 / dt_init); // how many timesteps in 1min, which is the actual simulation timestep
    double dx = 1.0;// / space_grid_controller; // space step in x direction, double to be consistent with other types
    double dy = 1.0;// / space_grid_controller; // space step in y direction

    // reaction rate
    double k_reac = 0.001;//1.0;//0.00001;//0.1;//0.105;//0.03;//.205; // reaction term

   bool cellsmove = false; // cells undergo active movement(e.g. chemotaxis, random movement)
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
    double speed_l = 0.14;// 2.5/5.0;//  // speed of a leader cell default: 0.14;, before I checked 0.8mu m/min, since 1 min is 0.05, 0.01 is 1/5 min
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



    // int n_seed = 0;
    // double diff_conc = 0.05;


    // domain growth parameters


    /*
    * strain rate
    * */






    VectorXd strain = VectorXd::Zero(length_x);




// cells hinder the growth


    /*
    * strain rate
    * */
    // G3 baseline
    double initial_strain = 0.01;// trying now 0.005, chemical ablation, reduction 0.005

    //piecewise constant, two parts
    // 1 part n_faster times faster than the other part
    double n_faster = 2.0;

    double thetasmall = 1.0; // first thetasmall is growing
    int theta1 = int(thetasmall * length_x);

    double alpha1;
    double alpha2;


    /*
     * Uncomment for cell-hindered growth (G4)
//     * */
    //check whether first part grows
//    double thetasmalltemp = thetasmall;
//    if (first_part_grows ==true){
//        double xvar = final_length/ (n_faster * double(length_x)*thetasmalltemp + double(length_x)*(1-thetasmalltemp)); // solve: 2 *xvar * length_x * thetasmall + xvar * length(1-thetasmall) = final_length
//
//        double ratio1 = n_faster * double(length_x)*thetasmalltemp * xvar / (double(length_x)*thetasmalltemp);
//
//        alpha1 = log (ratio1)/final_time;
//
//        double ratio2 =  double(length_x)*(1- thetasmalltemp) * xvar / (double(length_x)*(1- thetasmalltemp));
//
//        alpha2 = log (ratio2)/final_time;
//
//    }
//
//    else if(first_part_grows == false){
//        double xvar = final_length/ (double(length_x)*thetasmall + n_faster*double(length_x)*(1-thetasmall)); // solve: 2 *xvar * length_x * thetasmall + x * length(1-thetasmall) = final_length
//
//        double ratio1 = double(length_x)*thetasmall * xvar / (double(length_x)*thetasmall);
//
//        alpha1 = log (ratio1)/final_time;
//
//        double ratio2 =  double(length_x)*(1- thetasmall) * xvar*n_faster / (double(length_x)*(1- thetasmall));
//
//        alpha2 = log (ratio2)/final_time;
//    }
//
//
//
//    double strainInitial = alpha1 +domain_growth_par; // for all before I had +0.01
//
//    // first part it is linear growth
//    for (int i = 0; i < theta1; i++) {
//        strain(i) = alpha1+domain_growth_par;//linear_par * double(theta1) /
//        // double(space_grid_controller);//linear_par * (double(i) / double(space_grid_controller));
//    }
//
//    for (int i = theta1; i < length_x; i++) {
//        strain(i) = alpha2+domain_growth_par;// 0.002;//0.5 * linear_par * double(theta1) / double(space_grid_controller); // constant to where it was
//        //strain(i,j) = linear_par*theta1/(theta1- (theta2-1))*(i-(theta2-1)); // linearly decreasing, if I do this do not forget to change Gamma
//    }

//// end of G4


    // growth function

    // I will mainly use its first derivative with respect to x

    VectorXd Gamma_x = VectorXd::Zero(length_x);
    VectorXd Gamma_x_old = VectorXd::Zero(length_x);
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


//    // save data to plot chemoattractant concentration
//    ofstream outputinit("matrix_non_uniform0.csv");
//
//    outputinit << "x, y, z, u" << "\n" << endl;
//
//
//    for (int i = 0; i < length_x * length_y; i++) {
//        for (int j = 0; j < 4; j++) {
//            outputinit << chemo_3col(i, j) << ", ";
//        }
//        outputinit << "\n" << endl;
//    }


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
    ABORIA_VARIABLE(scaling, int, "scaling ") // stores the value that the cell position is scaled down to initial coordinates
    // the chain and then increasing integer if it is. If it is attached to a leader, it is 1, and then increasing order.
    // stores the distance to the closest neighbour, if less than thresold
    typedef Particles<std::tuple<radius, type, attached_to_id, direction, chain, chain_type, persistence_extent, same_dir_step,scaling>, 2> particle_type; // 2 stands for dimension

    // will use stored value of the position of a particle
    typedef particle_type::position position;

    // initialise the number of particles
    particle_type particles(N*6); // times 6 beforehand

    // initialise random number generator for particles entering the domain, appearing at the start in x and uniformly in y
    std::default_random_engine gen;
    std::uniform_real_distribution<double> uniform(cell_radius, length_y - 1 - cell_radius);


    /*
     * compact initialisation of particles
     */

    for (int i = 0; i < N; ++i) {


        get<radius>(particles[i]) = cell_radius;
        get<type>(particles[i]) = 0; // initially all cells are leaders

//                get<position>(particles[i]) = vdouble2(cell_radius, (i + 1) * double(length_y - 1) / double(N) -
//                                                            0.5 * double(length_y - 1) /
//                                                            double(N)); // x=2, uniformly in y


        // equidistributed cells
        //IMPORTANT UNCOMMENT
        get<position>(particles[i+5*N]) = vdouble2(cell_radius, (i + 1) * double(length_y - 1) / double(N) -
                                                            0.5 * double(length_y - 1) /
                                                            double(N)); // x=2, uniformly in y

        get<position>(particles[i+4*N]) = vdouble2(cell_radius + 50, (i + 1) * double(length_y - 1) / double(N) -
                                                            0.5 * double(length_y - 1) /
                                                            double(N)); // x=2, uniformly in y

        get<position>(particles[i+ 3*N]) = vdouble2(cell_radius + 100, (i + 1) * double(length_y - 1) / double(N) -
                                                            0.5 * double(length_y - 1) /
                                                            double(N)); // x=2, uniformly in y

        get<position>(particles[i + 2*N]) = vdouble2(cell_radius + 150, (i + 1) * double(length_y - 1) / double(N) -
                                                            0.5 * double(length_y - 1) /
                                                            double(N)); // x=2, uniformly in y

        get<position>(particles[i + N]) = vdouble2(cell_radius + 200, (i + 1) * double(length_y - 1) / double(N) -
                                                            0.5 * double(length_y - 1) /
                                                            double(N)); // x=2, uniformly in y

        get<position>(particles[i]) = vdouble2(cell_radius + 250, (i + 1) * double(length_y - 1) / double(N) -
                                                            0.5 * double(length_y - 1) /
                                                            double(N)); // x=2, uniformly in y

        get<persistence_extent>(particles[i]) = 0;
        get<same_dir_step>(particles[i]) = 0;

//        // for fixed cells
//        if (i==0){
//            get<position>(particles[i]) = vdouble2(cell_radius, (i + 1) * double(length_y - 1) / double(N) -
//                                                                0.5 * double(length_y - 1) /
//                                                                double(N)); //
//        }
//        if(i==1) {
//            get<position>(particles[i]) = vdouble2(170, (i + 1) * double(length_y - 1) / double(N) -
//                                                        0.5 * double(length_y - 1) /
//                                                        double(N)); //
//        }
//
//        if(i==2){
//            get<position>(particles[i]) = vdouble2(339, (i + 1) * double(length_y - 1) / double(N) -
//                                                        0.5 * double(length_y - 1) /
//                                                        double(N)); //
//        }




    }

    // initialise neighbourhood search, note that the domain will grow in x direction, so I initialise larger domain
    particles.init_neighbour_search(vdouble2(0, 0), 5 * vdouble2(length_x, length_y), vbool2(false, false));

    // save particles before they move

    //  vtkWriteGrid("particles", t, particles.get_grid(true));

    // initialise random number generator to obtain random number between 0 and 2*pi
    std::default_random_engine gen1;
    gen1.seed(t * n_seed); // choose different seeds to obtain different random numbers
    std::uniform_real_distribution<double> uniformpi(0, 2 * M_PI);


 //        ofstream outputtrackL("speed1p5track_leadTheta" + to_string(thetasmall) + "FINALnseed" + to_string(n_seed) + ".csv");
////        ofstream outputtrackL2("track_lead2theta" + to_string(thetasmall) + "nseed"+ to_string(n_seed) + ".csv");
////        ofstream outputtrackL3("track_lead3theta" + to_string(thetasmall) + "nseed"+ to_string(n_seed) + ".csv");
////        ofstream outputtrackL4("track_lead4theta" + to_string(thetasmall) + "nseed"+ to_string(n_seed) + ".csv");
////        ofstream outputtrackL5("track_lead5theta" + to_string(thetasmall) + "nseed"+ to_string(n_seed) + ".csv");
//
////    ofstream outputtrackF("CORRECTV2FOLnongrowingnseed"+ to_string(n_seed) + ".csv");
//          ofstream outputtrackF("speed1p5track_folTheta" + to_string(thetasmall) + "FINALnseed"+ to_string(n_seed) + ".csv");
////        ofstream outputtrackF2("track_fol2theta" + to_string(thetasmall) + "nseed"+ to_string(n_seed) + ".csv");
////        ofstream outputtrackF3("track_fol3theta" + to_string(thetasmall) + "nseed"+ to_string(n_seed) + ".csv");
////        ofstream outputtrackF4("track_fol4theta" + to_string(thetasmall) + "nseed"+ to_string(n_seed) + ".csv");
////        ofstream outputtrackF5("track_fol5theta" + to_string(thetasmall) + "nseed"+ to_string(n_seed) + ".csv");


    //for each timestep
    while (t < final_time) {



//              insert new cells
//

        //if (particles.size() < 50){

        bool free_position = false;
        particle_type::value_type f;
        //get<radius>(f) = cell_radius;


        get<position>(f) = vdouble2(cell_radius, uniform(gen)); // x=2, uniformly in y
        free_position = true;
        /*
         * loop over all neighbouring leaders within "dem_diameter" distance
         */
        for (auto tpl = euclidean_search(particles.get_query(), get<position>(f), diameter); tpl != false; ++tpl) {

            vdouble2 diffx = tpl.dx();

            if (diffx.norm() < diameter) {
                free_position = false;
                break;
            }
        }

        // our assumption that all new cells are followers
        get<type>(f) = 1;

        // this is if I only have a certain number of cells
        //if (t < 2) { // this

        if (free_position) {
            get<chain>(f) = 0;
            get<chain_type>(f) = -1;
            get<attached_to_id>(f) = -1;
            particles.push_back(f);
        }
       // } // this, less than 50 cells end
        particles.update_positions();






        /*
         * Remove some cells half time
         * */

        //cout << " t" << t << endl;
//        if (t >27.01 && t < 27.03){
//            cout << "in here " << endl;
//            for (auto p : particles) {
//                if (get<position>(p)[0] < 200.0 && get<position>(p)[0] > 150.000) { //before was 350 and 250
//                    get<alive>(p) = false;
//                }
//            }
//            particles.update_positions();
//        }




        t = t + dt;

        counter = counter + 1;

//


        /*
 *
 * Domain growth, and update cell positions
 */

        // CELL-INDUCED GROWTH, UPDATE STRAIN RATE BASED ON CELL POSITIONS

        double leaders_pos = 0.0;
        double av_lead = 0.0;
        vdouble2 lead_coord;
//        // position of leaders, take the average of leaders
//        for (int i = 0; i < N; i++) {
//            //cout << "position leaders " << get<position>(particles[i]) << endl;
//            lead_coord = get<position>(particles[i]);
//            leaders_pos = leaders_pos + lead_coord[0];
//        }
//        av_lead = leaders_pos / double(N); // not scaled to initial coordinates

//         if I only take the cells that is the furthest ahead
        for (int i = 0; i < N; i++) {
            lead_coord = get<position>(particles[i]);
            if (lead_coord[0] > leaders_pos) {
                av_lead = lead_coord[0];
                leaders_pos = lead_coord[0];
            }
        }

        int j = 0;
        while (av_lead  - 100 > Gamma_old(j)) {// // av_lead- 100.0> Gamma_old(j) (when with delay, G2) av_lead + 100.0 when ahead
            value = j;
            j = j + 1;

        }

        theta1 = value;

//    // linearly decreasing (G3) G4 if strain(i) non-zero for i>theta1


//        if (t < 40.5){
//            for (int i = 0; i < theta1 + 50; i++) {
//                strain(i) = domain_growth_par * i / theta1 + initial_strain; //when there is baseline growth // 0.035*(theta1-i)/theta1; // linearly increasing to the NT
//            }
//            for (int i = theta1 +50; i < length_x-1; i++) {
//                strain(i) = initial_strain; //when there is baseline growth // linearly increasing to the NT
//            }
//        }else{
        for (int i = 0; i < theta1; i++) {
            strain(i) = domain_growth_par +
                        initial_strain;// * (theta1-i) / theta1 + initial_strain; //when there is baseline growth // 0.035*(theta1-i)/theta1; // linearly increasing to the NT
        }
        for (int i = theta1; i < length_x - 1; i++) {
            strain(i) = initial_strain; //when there is baseline growth // linearly increasing to the NT
        }
//        }


        // Density dependent growth rate

        //split the domain into sections of length 100, the whole domain or from 0 to that length
//
//    const int domain_partition = int(Gamma(length_x - 1) / double(55)); // number of intervals of 50 \mu m
//
//
//    VectorXi parts = VectorXi::Zero(domain_partition); // integer with number of cells in particular part that plus one is proportions that are not in chains
//
//    double one_part = Gamma(length_x - 1) / double(domain_partition);
//
//
//    for (int i = 0; i < domain_partition; i++) {
//
//        for (int j = 0; j < particles.size(); j++) {
//            vdouble2 x = get<position>(particles[j]);
//            if (i * one_part < x[0] && x[0] < (i + 1) * one_part) {
//                parts(i) += 1;
//            }
//        }
//
//    }

//    cout << "one part " << one_part << endl;
//    cout << "parts " << parts << endl;
//    cout << "domain partition " << domain_partition << endl;
//    int k = 0;
//    int k_old = 0;
//    for (int i = 0; i < domain_partition; i++) {
//        if (parts(i) < 10){ //if less than ten cells in a section, then strain rate only as an initial strain rate
//                while ((i+1) * one_part > Gamma_old(k)) {// // av_lead-100.0> Gamma_old(j) (when with delay, G2)
//                    value = k;
//                    k = k + 1;
//                    }
//                for (int j=k_old;j<k;j++){
//                    strain(j) = initial_strain;// * (theta1-i) / theta1 + initial_strain; //when there is baseline growth // 0.035*(theta1-i)/theta1; // linearly increasing to the NT
//                }
//                k_old = k;
//            }
//
//
//        if (parts(i) < 20 && parts(i) >= 10){//if more than 10 but less than 20 cells in a section, then strain rate an initial strain rate +domain_growhth_par
//                while ((i+1) * one_part > Gamma_old(k)) {// // av_lead-100.0> Gamma_old(j) (when with delay, G2)
//                    value = k;
//                    k = k + 1;
//                    }
//                for (int j=k_old;j<k;j++){
//                    strain(j) = initial_strain+ domain_growth_par;// * (theta1-i) / theta1 + initial_strain; //when there is baseline growth // 0.035*(theta1-i)/theta1; // linearly increasing to the NT
//                }
//                k_old = k;
//        }
//        if ( parts(i) >= 20){ //if more than 20 cells in a section, then strain rate an initial strain rate + 2 * domain_growhth_par
//                while ((i+1) * one_part > Gamma_old(k)) {// // av_lead-100.0> Gamma_old(j) (when with delay, G2)
//                    value = k;
//                    k = k + 1;
//                    }
//                for (int j=k_old;j<k;j++){
//                    strain(j) = initial_strain+ 2* domain_growth_par;// * (theta1-i) / theta1 + initial_strain; //when there is baseline growth // 0.035*(theta1-i)/theta1; // linearly increasing to the NT
//                }
//                k_old = k;
//        }
//    }


        // end of density dependent domain growth


//    cout << "new strain " << endl;
//    cout << "strain rate " << strain << endl;



        // Cell-hindered (G5)
//        for (int i = 0; i < theta1; i++) {
//            strain(i) = strainInitial * (i) / theta1; // linearly increasing
//        }

//         G1 and G2
//        for (int i = 0; i < theta1; i++){
//            strain(i) = 0.025; // for all the examples before 0.025, when I had heaviside function
//        }

        if (counter == 1){
            for (int i = 0; i < length_x; i++) {
                Gamma_x_old(i) = exp(t * strain(i));
            }
        }
//
//        // update the strain rate
//        for (int i = 0; i < length_x; i++) {
//            Gamma_x(i) = Gamma_x_old(i) + exp((dt) * strain(i));
//        }

        for (int i = 0; i < length_x; i++) {
            Gamma_x(i) = Gamma_x_old(i)*exp(dt* strain(i));
        }

        Gamma_x_old = Gamma_x;

        // domain length


        /*
         * this is important and it will change based on strain rates, now since there is linear growth in the first section,
         * the first factor appears due to integration
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


//        for (int i = 0; i < theta1; i++) {
////            Gamma(i) =
////                    Gamma_initial * (double(i) / double(space_grid_controller)) * (Gamma_x(i)); // linearly increasing
//            Gamma(i) = (double(i)) * (Gamma_x(i));
//            cout << "old Gamma(i) " << i << " value " << Gamma(i) << endl;
//
//        }
//
//
//
//        // second part is constant
//        for (int i = theta1; i < length_x; i++) {
////
////            Gamma(i) = Gamma(theta1 - 1) +
////                       Gamma_initial * (double(i) / double(space_grid_controller) - double(theta1 - 1) /
////                                                                                    double(space_grid_controller)) *
////                       Gamma_x(i);
//
//            Gamma(i) = Gamma(theta1 - 1) + (double(i) - double(theta1 - 1)) * Gamma_x(i);
////
//            cout << "old Gamma(i) " << i << " value " << Gamma(i) << endl;
//            //Gamma(i,j) = ; // linearly decreasing, if I do this do not forget to change Gamma
//
//        }


        /*
         * arbitrary Gamma function
         * */
// old

        Gamma(0) = 0; // this is assumption, since I cannot calculate it
        //cout << "Gamma(0) " << 0 << " value " << Gamma(0) << endl;

        for (int i = 1; i < length_x; i++) {

            Gamma(i) = Gamma_x(i) * dx + Gamma(i - 1);

            //   cout << "Gamma(i) " << i << " value " << Gamma(i) << endl;
        }


//new

//        Gamma_temp(0) = 0; // this is assumption, since I cannot calculate it
//
//        double total_increase = 0;
//        double smallest_increase = 0;
//        double total_sum = Gamma_x(0);
//        double factor = 0; // sum all gamma over this factor is equal to total change
//
////            cout << "theta1 " << theta1 << endl;
////            cout << "Gamma_base " << Gamma_base(theta1) << endl;
//
//            // old version which worked
//            //   this is for general case
//            for (int i = 1; i < length_x; i++) {
//
//                Gamma_temp(i) = Gamma_x(i) * dx + Gamma_temp(i - 1);
//                total_sum = total_sum + Gamma_x(i);
//            }
//            total_increase = Gamma_temp(length_x - 1) - Gamma(length_x - 1);
//
//            factor = Gamma_temp(length_x-1) / total_increase;
//
//            for (int i = 1; i < length_x; i++) {
//                Gamma(i) = Gamma_old(i) + Gamma_temp(i) / factor;
//            }






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


        //physical aplabtion of cells, cells only start moving after 2 hours (t > 6), 4hrs (t>12)

        //if (t > 12) {
        /// update positions uniformly based on the domain growth

        vdouble2 x; // use variable x for the position of cells
        double x0 = 0;
        int pos;

        for (int i = 0; i < particles.size(); i++) {

            x = get<position>(particles[i]);

            // since I do not know how to do it for general case, I will do it for my specific

            // analytical for two different regions
            //            if (x[0] > Gamma_old(theta1 - 1)) {
            //
            //                x0 = (x[0] - Gamma_old(theta1 - 1)) *
            //                        ((Gamma(length_x - 1) - Gamma(theta1 - 1)) /
            //                     (Gamma_old(length_x - 1) - Gamma_old(theta1 - 1))) + Gamma(theta1 - 1);
            //
            //                get<position>(particles)[i] = vdouble2(x0,x[1]);
            //
            //            } else {
            //                get<position>(particles)[i] *= vdouble2((Gamma(theta1 - 1)) / (Gamma_old(theta1 - 1)),
            //                                                        1); // update position based on changes in Gamma
            //            }
            int j = 0;
            while (x[0] > Gamma_old(j)) {
                value = j;
                j = j + 1;
                //cout << "value " << value << endl;

            }


            get<scaling>(particles)[i] = value;

            //            x[0] = x[0] + Gamma(value)-Gamma_old(value);

            get<position>(particles)[i] += vdouble2(Gamma(value) - Gamma_old(value), 0);


        }
//
        // } // physical ablation
        Gamma_old = Gamma;


        /*
         * update chamo concentration
         *
         * */


        // internalisation rate



        // internalisation
        for (int i = 0; i < length_x; i++) {
            for (int jGamma_x_old = 0; j < length_y; j++) {
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
                                        ((1.0 / Gamma_x(i) + 1.0 / Gamma_x(i + 1)) * (chemo(i + 1, j) - chemo(i, j)) -
                                         (chemo(i, j) - chemo(i - 1, j)) * (1.0 / Gamma_x(i) + 1.0 / Gamma_x(i - 1))) +
                                        D * (chemo(i, j + 1) - 2 * chemo(i, j) + chemo(i, j - 1)) / (dy * dy) -
                                        (chemo(i, j) * lam / (2 * M_PI * cell_radius * cell_radius)) * intern(i, j) +
                                        chemo(i, j) * k_reac * (1 - chemo(i, j)) - strain(i) * chemo(i, j)) +
                                  chemo(i, j);
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
//




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

        // physical ablation of NC cells, 1hrs 6, 4hrs 12
        // if (t>12) {
        // update the position of all particles in a random order created above

//        for (int i=0; i < particles.size(); i++){
//            cout << get<position>(particles[i]) << endl;
//        }
        if (cellsmove == true) {

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

                // if two regions
                //                if (x[0] < Gamma(theta1 - 1)) {
                //                    x_in = x[0] * (theta1) / Gamma(theta1 - 1);
                //                    l_filo_x = l_filo_x_in * theta1 / Gamma(theta1 - 1);
                //                } else {
                //                    x_in = (x[0]-Gamma(theta1-1))* (length_x-theta1)/(Gamma(length_x-1)-Gamma(theta1-1)) + theta1;
                //                    l_filo_x = l_filo_x_in * (length_x - theta1) / (Gamma(length_x-1 ) - Gamma(theta1 - 1));
                //                }

                // for any Gamma

                x_in = get<scaling>(particles)[particle_id(j)];
                l_filo_x = l_filo_x_in * get<scaling>(particles)[particle_id(j)] /
                           Gamma(get<scaling>(particles)[particle_id(j)]);


                // if it is still in the process of moving in the same direction
                if (get<persistence_extent>(particles[particle_id(j)]) == 1) {


                    // two regions
                    //                    x += get<direction>(particles[particle_id(j)]);
                    //
                    //                    if (x[0] < Gamma(theta1 - 1)) {
                    //                        x_in = x[0] * (theta1) / Gamma(theta1 - 1);
                    //                        l_filo_x = l_filo_x_in * theta1 / Gamma(theta1 - 1);
                    //                    } else {
                    //                        x_in = (x[0]-Gamma(theta1-1))* (length_x-theta1)/(Gamma(length_x-1)-Gamma(theta1-1)) + theta1;
                    //                        l_filo_x = l_filo_x_in * (length_x - theta1) / (Gamma(length_x-1 ) - Gamma(theta1 - 1));
                    //                    }



                    bool free_position = true; // check if the neighbouring position is free

                    // check if there are other particles in the position where the particle wants to move
                    for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {
                        if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            free_position = false;
                        }
                    }

                    // check that the position they want to move to is free and not out of bounds
                    if (free_position && x[0] > cell_radius && x[0] < Gamma(length_x - 1) - cell_radius &&
                        (x[1]) > cell_radius &&
                        (x[1]) < length_y - 1 - cell_radius) {
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

                        // two regions
                        //                        if (x[0] < Gamma(theta1 - 1)) {
                        //                            x_in = x[0] * (theta1) / Gamma(theta1 - 1);
                        //                            l_filo_x = l_filo_x_in * theta1 / Gamma(theta1 - 1);
                        //                        } else {
                        //                            x_in = (x[0]-Gamma(theta1-1))* (length_x-theta1)/(Gamma(length_x-1)-Gamma(theta1-1)) + theta1;
                        //                            l_filo_x = l_filo_x_in * (length_x - theta1) / (Gamma(length_x-1 ) - Gamma(theta1 - 1));
                        //                        }

                        // for any Gamma

                        bool free_position = true; // check if the neighbouring position is free

                        // check if the position the particle wants to move is free
                        for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                            if (get<id>(*k) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                free_position = false;
                            }
                        }


                        // if the position they want to move to is free and not out of bounds, move that direction
                        if (free_position && x[0] > cell_radius && x[0] < Gamma(length_x - 1) - cell_radius &&
                            (x[1]) > cell_radius &&
                            (x[1]) < length_y - 1 - cell_radius) {
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

                        // for specific case
                        //                        if (x[0] < Gamma(theta1 - 1)) {
                        //                            x_in = x[0] * (theta1) / Gamma(theta1 - 1);
                        //                            l_filo_x = l_filo_x_in * theta1 / Gamma(theta1 - 1);
                        //                        } else {
                        //                            x_in = (x[0]-Gamma(theta1-1))* (length_x-theta1)/(Gamma(length_x-1)-Gamma(theta1-1)) + theta1;
                        //                            l_filo_x = l_filo_x_in * (length_x - theta1) / (Gamma(length_x-1 ) - Gamma(theta1 - 1));
                        //                        }



                        bool free_position = true; // check if the neighbouring position is free

                        // if this loop is entered, it means that there is another cell where I want to move
                        for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                            if (get<id>(*k) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                free_position = false;
                            }
                        }


                        // update the position if the place they want to move to is free and not out of bounds
                        if (free_position && x[0] > cell_radius && x[0] < Gamma(length_x - 1) - cell_radius &&
                            (x[1]) > cell_radius &&
                            (x[1]) < length_y - 1 - cell_radius) {
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

                // fro two regions
                //                if (x[0] < Gamma(theta1 - 1)) {
                //                    x_in = x[0] * (theta1) / Gamma(theta1 - 1);
                //                    l_filo_x = l_filo_x_in * theta1 / Gamma(theta1 - 1);
                //                } else {
                //                    x_in = (x[0]-Gamma(theta1-1))* (length_x-theta1)/(Gamma(length_x-1)-Gamma(theta1-1)) + theta1;
                //                    l_filo_x = l_filo_x_in * (length_x - theta1) / (Gamma(length_x-1 ) - Gamma(theta1 - 1));
                //                }

                //for any Gamma

                x_in = get<scaling>(particles)[j];
                l_filo_x = l_filo_x_in * get<scaling>(particles)[particle_id(j)] /
                           Gamma(get<scaling>(particles)[particle_id(j)]);

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

                    // for two regions
                    //                    if (x_chain[0] < Gamma(theta1 - 1)) {
                    //                        x_in_chain = x_chain[0] * (theta1) / Gamma(theta1 - 1);
                    //                        l_filo_x = l_filo_x_in * theta1 / Gamma(theta1 - 1);
                    //                    } else {
                    //                        x_in_chain = (x_chain[0]-Gamma(theta1-1))* (length_x-theta1)/(Gamma(length_x-1)-Gamma(theta1-1)) + theta1;
                    //
                    //                        l_filo_x = l_filo_x_in * (length_x - theta1) / (Gamma(length_x-1 ) - Gamma(theta1 - 1));
                    //                    }

                    //for any Gamma


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
                    if (free_position && x_chain[0] > cell_radius && x_chain[0] < Gamma(length_x - 1) - cell_radius &&
                        (x_chain[1]) > cell_radius &&
                        (x_chain[1]) < length_y - 1 - cell_radius) {
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
                    // l_filo_x not scaled because I am looking at cell positions which are not scaled
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


                        // for two regions
                        //                        if (x_chain[0] < Gamma(theta1 - 1)) {
                        //                            x_in_chain = x_chain[0] * (theta1) / Gamma(theta1 - 1);
                        //                        } else {
                        //                            x_in_chain = (x_chain[0]-Gamma(theta1-1))* (length_x-theta1)/(Gamma(length_x-1)-Gamma(theta1-1)) + theta1;
                        //                            l_filo_x = l_filo_x_in * (length_x - theta1) / (Gamma(length_x-1 ) - Gamma(theta1 - 1));
                        //                        }

                        //for any Gamma

                        //x_in_chain = get<scaling>(particles)[particle_id(j)];



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
                            x_chain[0] < Gamma(length_x - 1) - cell_radius && (x_chain[1]) > cell_radius &&
                            (x_chain[1]) < length_y - 1 - cell_radius) {
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

                        //                        while (((x_in + sin(random_angle) * l_filo_x_in) < 0 ||
                        //                                ((x_in + sin(random_angle) * l_filo_x_in)) >
                        //                                length_x - 1 || (x[1] + cos(random_angle) * l_filo_y) < 0 ||
                        //                                (x[1] + cos(random_angle) * l_filo_y) > length_y - 1)) {
                        //                            random_angle = uniformpi(gen1);
                        //                        }

                        x += speed_f * vdouble2(sin(random_angle), cos(random_angle));


                        // FOR TWO REGIONS
                        //                        if (x[0] < Gamma(theta1 - 1)) {
                        //                            x_in = x[0] * (theta1) / Gamma(theta1 - 1);
                        //                        } else {
                        //                            x_in = (x[0]-Gamma(theta1-1))* (length_x-theta1)/(Gamma(length_x-1)-Gamma(theta1-1)) + theta1;
                        //                        }

                        //for any Gamma

                        // x_in = get<scaling>(particles)[j];



                        bool free_position = true; // check if the neighbouring position is free

                        // check if the position the cells want to move to is free
                        for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                            if (get<id>(*k) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                free_position = false;
                            }
                        }

                        // if the position they want to move to is free and not out of bounds, move to that position
                        if (free_position && x[0] > cell_radius && x[0] < Gamma(length_x - 1) - cell_radius &&
                            (x[1]) > cell_radius &&
                            (x[1]) < length_y - 1 - cell_radius) {
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
        // } // physical ablation
    } //cells move
        // update positions
        particles.update_positions();


        if (counter % 100 == 0) {

//
//            #ifdef HAVE_VTK
//                        vtkWriteGrid("ToSaveOLDCellupdateminus100domgrowthpar0p01Cell_inducedCELLS", t, particles.get_grid(true));
//            #endif
//            //
//            //
//            ////
//           //  //ofstream output("matrix_FIRST_025theta" + to_string(int(round(t))) + ".csv");
//            ofstream output("ToSaveOLDcellupdateminus100domgrowthpar0p01Cell_inducedMATRIX" + to_string(int(t)) + ".csv");
//

#ifdef HAVE_VTK
            vtkWriteGrid("EQuidistrdomgrowthrate0p018CELLS", t, particles.get_grid(true));
#endif
            //
            //
            ////
            //  //ofstream output("matrix_FIRST_025theta" + to_string(int(round(t))) + ".csv");
            ofstream output("EQuidistrdomgrowthrate0p018Matrix" + to_string(int(t)) + ".csv");




            output << "x, y, z, u" << "\n" << endl;



            for (int i = 0; i < length_x * length_y; i++) {
                for (int j = 0; j < 4; j++) {
                    output << chemo_3col(i, j) << ", ";
                }
                output << "\n" << endl;
            }


        }
    }

//    /*
// * return the density of cells in domain_partition parts of the domain
// */
    const int domain_partition = int(Gamma(length_x - 1) / double(55));; // number of intervals of 50 \mu m


    VectorXi proportions = VectorXi::Zero(domain_partition); // integer with number of cells in particular part that plus one is proportions that are not in chains

    double one_part = Gamma(length_x - 1) / double(domain_partition);


    for (int i = 0; i < domain_partition; i++) {

        for (int j = 0; j < particles.size(); j++) {
            vdouble2 x = get<position>(particles[j]);
            if (i * one_part < x[0] && x[0] < (i + 1) * one_part) {
                proportions(i) += 1;
            }
        }

    }

    //print the proportion that break
    double pro_break = 0.0; // the stream did not break

    int followers_not_in_chain = 0;

    for (int i = 0; i < particles.size(); ++i){
        if (get<chain>(particles[i]) == 0){
            if (get<type>(particles[i]) == 1){
                followers_not_in_chain += 1; // add to coung
            }
        }
    }

    pro_break = double(followers_not_in_chain)/(double(particles.size()-N));

    //position of leaders
    for (int i = 0; i<N; i++){
        cout << "position leaders " << get<position>(particles[i]) << endl;
    }
    cout << "domain growth par = " << domain_growth_par << endl;

    cout << "domain length " <<  Gamma(length_x-1) << endl;
    cout << "prop break " << pro_break << endl;

    // Gamma
//    cout << "Gamma(length_x-1) " << Gamma(length_x-1) << endl;




    return proportions;

}


/*
 * main for proportions in different sections
 */


// parameter analysis
int main() {

    const int number_parameters = 1; // parameter range
    const int sim_num = 1;

double domain_growth_par;
    double speed;
    // n would correspond to different seeds
    // parallel programming
    #pragma omp parallel for
        for (int k = 0; k < 1; k++){
            //domain_growth_par = 0.006 + 0.002*double(k);
            domain_growth_par = 0.018;//0.018;//0.015 for cell hindered, speed 0.14*1.5!!!! //cellinduced: 0.018; // cellhindered 0.006; //cellinduced: 0.018; //for G3: 0.01 + 0.0013*double(k);, for G4: 0.006 + 0.0005 * double(k)  // test 0.018
            for (int n = 0; n < sim_num; n++) {
                proportions(0.01, n, domain_growth_par); //0.01 is the control
            }

        }


}
