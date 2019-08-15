//
// Created by giniunaite on 06/08/19.
//

//
// Created by giniunaite on 14/05/19.
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




VectorXd proportions(double diff_conc, int n_seed, double domain_growth_par, double speed) {
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

    double first_profile  = 60.0;// growth profile changes at this time
// parameters for the dynamics of chemoattractant concentration

    double D = 2.0;//0.00001;//0.05; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0.0; // initialise time
    double dt = 0.1; // time step
    double dt_init = dt;
    int number_time = int(1 / dt_init); // how many timesteps in 1min, which is the actual simulation timestep
    double dx = 1.0;// / space_grid_controller; // space step in x direction, double to be consistent with other types
    double dy = 1.0;// / space_grid_controller; // space step in y direction

    // reaction rate
    double k_reac = 1.0;//0.00001;//0.1;//0.105;//0.03;//.205; // reaction term


    // cell parameters

    double cell_radius = 7.5;//0.5; // radius of a cell
    const double diameter =
            2 * cell_radius; // diameter of a cell
    const size_t N = 1; // initial number of cells
    double l_filo_y = 27.5;//2; // sensing radius, filopodia + cell radius
    double l_filo_x = 27.5; // sensing radius, it will have to be rescaled when domain grows
    double l_filo_x_in = l_filo_x; // this value is used for rescaling when domain grows based on initial value
    double l_filo_max = 45; // this is the length when two cells which were previously in a chain become dettached
    //double diff_conc = 0.1; // sensing threshold, i.e. how much concentration has to be bigger, so that the cell moves in that direction
    int freq_growth = 1; // determines how frequently domain grows (actually not relevant because it will go every timestep)
    int insertion_freq = 1; // determines how frequently new cells are inserted, regulates the density of population
    double speed_l = 0.4*speed;//  speed of a leader cell (0.14, in this case since dt =0.1, 1.4 if with random movements). dt 2min/ Before I had 0.8
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


    //piecewise constant, two parts
    // 1 part n_faster times faster than the other part
    double n_faster = 2.0;

    double thetasmall = 1.0; // first thetasmall is growing
    int theta1 = int(thetasmall * length_x);

    double alpha1;
    double alpha2;



    // G3 baseline
    double initial_strain = 0.0; // 0.01 now

    /*
     * Uncomment for cell-hindered growth (G4)
     * */
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



    double strainInitial = alpha1 +domain_growth_par; // for all before I had +0.01

    // first part it is linear growth
    for (int i = 0; i < theta1; i++) {
        strain(i) = alpha1+domain_growth_par;//linear_par * double(theta1) /
        // double(space_grid_controller);//linear_par * (double(i) / double(space_grid_controller));
    }

    for (int i = theta1; i < length_x; i++) {
        strain(i) = alpha2+domain_growth_par;// 0.002;//0.5 * linear_par * double(theta1) / double(space_grid_controller); // constant to where it was
        //strain(i,j) = linear_par*theta1/(theta1- (theta2-1))*(i-(theta2-1)); // linearly decreasing, if I do this do not forget to change Gamma
    }

// end of G4

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
    particle_type particles(N);

    // initialise random number generator for particles entering the domain, appearing at the start in x and uniformly in y
    std::default_random_engine gen;
    std::uniform_real_distribution<double> uniform(cell_radius, length_y - 1 - cell_radius);


    /*
     * compact initialisation of particles
     */

    for (int i = 0; i < N; ++i) {


        get<radius>(particles[i]) = cell_radius;
        get<type>(particles[i]) = 0; // initially all cells are leaders

        //IMPORTANT UNCOMMENT
        get<position>(particles[i]) = vdouble2(cell_radius, (i + 1) * double(length_y - 1) / double(N) -
                                                            0.5 * double(length_y - 1) /
                                                            double(N)); // x=2, uniformly in y
        get<persistence_extent>(particles[i]) = 0;
        get<same_dir_step>(particles[i]) = 0;





    }

    // initialise neighbourhood search, note that the domain will grow in x direction, so I initialise larger domain
    particles.init_neighbour_search(vdouble2(0, 0), 5 * vdouble2(length_x, length_y), vbool2(false, false));

    // save particles before they move

    //  vtkWriteGrid("particles", t, particles.get_grid(true));

    // initialise random number generator to obtain random number between 0 and 2*pi
    std::default_random_engine gen1;
    gen1.seed(t * n_seed); // choose different seeds to obtain different random numbers
    std::uniform_real_distribution<double> uniformpi(0, 2 * M_PI);



    //for each timestep
    while (t < final_time){





        t = t + dt;

        counter = counter + 1;


        /*
 *
 * Domain growth, and update cell positions
 */

        // CELL-INDUCED GROWTH, UPDATE STRAIN RATE BASED ON CELL POSITIONS

        double leaders_pos = 0.0;
        double av_lead = 0.0;
        vdouble2 lead_coord;
        //position of leaders
        for (int i = 0; i<N; i++){
            //cout << "position leaders " << get<position>(particles[i]) << endl;
            lead_coord = get<position>(particles[i]);
            leaders_pos = leaders_pos + lead_coord[0];
        }
        av_lead = leaders_pos/double(N); // not scaled to initial coordinates

        int j = 0;
        while (av_lead  > Gamma_old(j)){// // av_lead-100.0> Gamma_old(j) (when with delay, G2)

            value = j;
            j = j+1;
            if (j > 341){
                break;
            }

        }

        theta1 = value;


        // when change uncomment and
        if (t > first_profile) {
        // linearly increasing rate (G3)

//        for (int i = 0; i < theta1; i++) {
//            strain(i) = domain_growth_par * (theta1 - i) / theta1 + initial_strain; //when there is baseline growth // linearly increasing to the NT
//        }
//        for (int i = theta1; i < length_x-1; i++) {
//                strain(i) = initial_strain; //when there is baseline growth // linearly increasing to the NT
//            }

        }

        if (t < first_profile){
//            // linearly increasing rate (G3)
//
            for (int i = 0; i < theta1; i++) {
                strain(i) = domain_growth_par * (theta1 - i) / theta1 + initial_strain; //when there is baseline growth // linearly increasing to the NT
            }
            for (int i = theta1; i < length_x-1; i++) {
                strain(i) = initial_strain; //when there is baseline growth // linearly increasing to the NT
            }

            // Cell-hindered (G4)
//            for (int i = 0; i < theta1; i++){
//                strain(i) = strainInitial* (i)/theta1; // linearly decreasing
//            }
        }


        // G1 and G2
//        for (int i = 0; i < theta1; i++){
//            strain(i) = 0.025; // for all the examples before 0.025, when I had heaviside function
//        }
//
//
        // update the strain rate
        for (int i = 0; i < length_x; i++) {
            Gamma_x(i) = exp(t * strain(i));
        }


        /*
         * arbitrary Gamma function
         * */


        Gamma(0) = 0; // this is assumption, since I cannot calculate it
        //cout << "Gamma(0) " << 0 << " value " << Gamma(0) << endl;

        if (t < first_profile ){
            for (int i = 1; i < length_x;i++){

                Gamma(i) = Gamma_x(i) * dx + Gamma(i-1);
            }
            Gamma_base = Gamma;
        }





        Gamma_temp(0) = 0; // this is assumption, since I cannot calculate it

        double total_increase = 0;
        double smallest_increase = 0;
        double total_sum = Gamma_x(0);
        double factor = 0; // sum all gamma over this factor is equal to total change


        if (t > first_profile ){


            // old version which worked
            //   this is for general case
            for (int i = 1; i < length_x;i++){

                Gamma_temp(i) = Gamma_x(i) * dx + Gamma_temp(i-1);
                total_sum = total_sum +Gamma_x(i);
            }
            total_increase = Gamma_temp(length_x-1) - Gamma(length_x-1);

            factor = total_sum/total_increase;


            for (int i = 1; i < length_x;i++){
                Gamma(i) = Gamma_old(i) + Gamma_temp(i)/factor;
            }

        }




        // I need Gamma_t for cos verification as well

        for (int i = 0; i < length_x; ++i) {

            Gamma_t(i) = (Gamma(i) - Gamma_old(i)) / dt;

        }



        /// update positions uniformly based on the domain growth

        vdouble2 x; // use variable x for the position of cells
        double x0=0;
        int pos;

        for (int i = 0; i < particles.size(); i++) {

            x = get<position>(particles[i]);


            int j = 0;
            while (x[0]> Gamma_old(j)){

                value = j;
                j = j+1;
                if (j > 341){
                    break;
                }
            }


            get<scaling>(particles)[i] = value;

            vdouble2 x; // use variable x for the position of cells
            x = get<position>(particles[0]);

            if (x[0] < Gamma(length_x-1) - (Gamma(value)-Gamma_old(value),0)){
                get<position>(particles)[i] += vdouble2(Gamma(value)-Gamma_old(value),0);
            }
        }


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




            // if a particle is a leader


            vdouble2 x; // use variable x for the position of cells
            x = get<position>(particles[particle_id(j)]);

            if (x[0] < Gamma(length_x-1) - speed_l ){ // if it is not too far ahead
                get<position>(particles)[particle_id(j)] += speed_l *
                                                            vdouble2(1, 0);; // update if nothing is in the next position


            } // update if nothing is in the next position


        }

        // update positions
        particles.update_positions();


//if ( counter % 10 == 0) {
//#ifdef HAVE_VTK
//    vtkWriteGrid("Cell", t, particles.get_grid(true));
//#endif
//
//
//    //
//    //ofstream output("matrix_FIRST_025theta" + to_string(int(round(t))) + ".csv");
//    ofstream output("MATRIX" + to_string(int(t)) + ".csv");
//
//
//    output << "x " << "\n" << endl;
//
//
//
//    //output << "x, y, z, u" << "\n" << endl;
//
//
//    for (int i = 0; i < length_x * length_y; i++) {
//        for (int j = 0; j < 4; j++) {
//            output << chemo_3col(i, j) << ", ";
//        }
//        output << "\n" << endl;
//    }
//
//}


    }



    //position of leaders
//    for (int i = 0; i<N; i++){
//        cout << "position leaders " << get<position>(particles[i]) << endl;
//    }
//    cout << "domain growth par = " << domain_growth_par << endl;
//
//    cout << "domain length " <<  Gamma(length_x-1) << endl;
//    cout << "speed " << speed << endl;

    // Gamma
//    cout << "Gamma(length_x-1) " << Gamma(length_x-1) << endl;

VectorXd proportions = VectorXd::Zero(4);
vdouble2 x =get<position>(particles[0]);

    proportions(0) = speed; // position
    proportions(1) = domain_growth_par; // position
    proportions(2) = x[0]; // position
    proportions(3) = Gamma(length_x-1); // position


    return proportions;

}


/*
 * main for proportions in different sections
 */


// parameter analysis
int main() {

    const int number_parameters = 1; // parameter range
    const int sim_num = 10;

    double domain_growth_par;
    double speed;
    int n =1;

    VectorXd prop = VectorXd::Zero(4);
    VectorXd numbers = VectorXd::Zero(4);

    // n would correspond to different seeds
    // parallel programming
    ofstream output3("CorrectCellInducedSMALLERSPEEDData.csv"); // relabel G4, G3

    for (int k = 0; k < 30; k++){
        domain_growth_par = 0.01 + 0.0013*double(k); //for G3: 0.01 + 0.0013*double(k);, for G4: 0.006 + 0.0005 * double(k)
        for (int j = 0; j < 30; j++){
            speed = 1.0 + 0.1*double(j);
                numbers =  proportions(0.05,n, domain_growth_par, speed); //just to know what the length is
            for (int count = 0; count < 4; count ++){
                output3 << numbers(count)  << ", ";
            }
            output3 << "\n" << endl;
        }

    }


}
