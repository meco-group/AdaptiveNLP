#include "core.hpp"
#include "makeBuildingBlocks.hpp"
#include "plotter.hpp"
#include "adaptiveCorridorExampleHelpers.hpp"
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

int main(){
    // if set to false, the MPC example is ran just once,
    // if set to true, the MPC example is ran for different values of N
    const bool COMPUTE_SCALING_RESULTS = false;

    //////////////////////////////////////
    // general definition of parameters //
    //////////////////////////////////////

    // Horizon length
    double T = 5.0;

    // Maximum number of time-steps AdaptiveNLP will take into account
    int Nmax = 200;

    // Maximum number of instances of an extra constraint
    int max_nb_extra_instances = 10;

    // Number of discretization intervals
    std::vector<int> nb_intervals_list;
    if (COMPUTE_SCALING_RESULTS){
        nb_intervals_list = {8, 15, 25, 35, 45};
    } else {
        nb_intervals_list = {8};
    }

    // Number of intervals discarded in the first MPC iteration (shorter
    // horizon that will be extended as the vehicle starts to move)
    int nb_reduced_horizon = 0;

    // Number of time-steps involved in every interval (nk)
    int nb_steps = 3;

    // Initial state
    std::vector<double> x0 = {-9.8, 1.0, 0.1, -0.1};

    // Viewing radius
    double view_radius = 10;

    // x-coordinates of obstacles
    std::vector<double> obstacle_x = 
        {40.0, 35.0, 45.0, 47.0, 46.0, 55.0, 50.5, 49};

    // y-coordinates of obstacles
    std::vector<double> obstacle_y = 
        {2.5,  1.2,  9.5,  7.0,  0.0,  8.5,  3.0,  2.5};

    // radii of obstacles
    std::vector<double> obstacle_r = 
        {2.0,  1.0,  3.3,  1.2,  1.9,  4.8,  0.8,  1.3};

    // x-coordinates of people
    std::vector<double> low_speed_pos_x = 
        {3.9, 14, 23};

    // y-coordinates of people
    std::vector<double> low_speed_pos_y = 
        {1.2, 3.0, 5.2};

    // radii of safety zones around people
    std::vector<double> low_speed_r = 
        {2.7, 2.7, 2.7};

    // allowed velocity at the position of a person
    double low_speed_min = 0.0;

    ///////////////////////////////
    // Prepare for multiple runs //
    ///////////////////////////////
    // number of runs to be performed
    int nb_runs = 1;

    // storage for latest timings results
    std::vector<std::vector<double>> tt;

    // storage for t_solve, t_update and t_total for every case
    std::vector<double> solve_timings_adaptive = {};
    std::vector<double> update_timings_adaptive = {};
    std::vector<double> total_timings_adaptive = {};
    std::vector<double> solve_timings_casadi_1 = {};
    std::vector<double> update_timings_casadi_1 = {};
    std::vector<double> total_timings_casadi_1 = {};
    std::vector<double> solve_timings_casadi_2 = {};
    std::vector<double> update_timings_casadi_2 = {};
    std::vector<double> total_timings_casadi_2 = {};

    // perform every run
    for (int run_counter = 0; run_counter < nb_runs; run_counter++){
    cout<<endl<<endl<<endl<<endl;
    cout<<"STARTING RUN "<<run_counter<<endl;
    cout<<endl<<endl<<endl<<endl;

    /////////////////////////////
    // create object instances //
    /////////////////////////////

    // create builing blocks for the specific NLP and extract the corridors
    std::vector<double> corridor_1;
    std::vector<double> corridor_2;
    BuildingBlocks blocks = makeBlocksCorridors(nb_steps, corridor_1, 
                                                corridor_2);
    int nx = blocks.get_nx();
    int nu = blocks.get_nu();

    // create plotter object used for writing output data to files
    Plotter plotter = Plotter(blocks.getFreeTime(), blocks.get_nx(), 
                              blocks.get_nu(), view_radius);
    std::string file_path = __FILE__;
    std::string dir_path = file_path.substr(0, file_path.rfind("/"));
    plotter.setOutputFolder(dir_path + "/plotting_data/");

    // add all obstacles, people and corridors to the plotter object
    for (int i = 0; i < obstacle_x.size(); i++){
        plotter.addCircle(obstacle_x[i], obstacle_y[i], obstacle_r[i]);
        plotter.setCircleVisibility(i, false);
    }
    for (int i = 0; i < low_speed_pos_x.size(); i++){
        plotter.addCircle(low_speed_pos_x[i], low_speed_pos_y[i], 
                          low_speed_r[i]);
        plotter.setCircleVisibility(obstacle_x.size() + i, false);
        plotter.setCircleColor(obstacle_x.size() + i, "b");
    }
    plotter.addCorridor(corridor_1);
    plotter.addCorridor(corridor_2);
    
    // create adaptive NLP object
    for (int nb_intervals : nb_intervals_list){
    double dt = T/((nb_intervals)*(nb_steps-1));
    double T_reduced = T - nb_reduced_horizon*(nb_steps-1)*dt;

    // number of MPC iterations
    int nb_iterations = int(40.625/dt);
    cout<<"nb_iterations = "<<nb_iterations<<endl;

    AdaptiveNLP adaptiveNLP = AdaptiveNLP(blocks, T_reduced, Nmax, 
                                          max_nb_extra_instances);

    // initialize some discretization intervals
    std::vector<int> nks(nb_intervals-nb_reduced_horizon, nb_steps);
    adaptiveNLP.initTimeSteps(nks);

    // add corridor 1 constraint to all time instances
    std::vector<int> all(adaptiveNLP.getN());
    for (int i = 0; i < all.size(); i++){
        all[i] = i + (i >= adaptiveNLP.getFinalInd());
    }
    adaptiveNLP.addExtraConstraint(all, {corridor_1});

    // create helper object used to change constraints
    AdaptiveCorridorHelper helper = AdaptiveCorridorHelper(adaptiveNLP, 
                                                           view_radius, 
                                                           low_speed_min);

    // perform adaptive approach
    tt = helper.performAdaptiveLoop(adaptiveNLP, plotter, x0, 
                               corridor_1, corridor_2, nx, nu, nb_steps, 
                               nb_iterations, dt, nb_intervals, obstacle_x, 
                               obstacle_y, obstacle_r, low_speed_pos_x, 
                               low_speed_pos_y, low_speed_r);

    // store timing results
    for (int i = 0; i < tt[0].size(); i++){
        solve_timings_adaptive.push_back(tt[0][i]);
        total_timings_adaptive.push_back(tt[1][i]);
        update_timings_adaptive.push_back(tt[1][i] - tt[0][i]);
    }

    // perform case CasADi Opti 1 (ref)
    tt = helper.performCasadiLoop(blocks, plotter, x0, corridor_1, 
                             corridor_2, nx, nu, nb_steps, nb_iterations, dt, 
                             T, nb_intervals, obstacle_x, obstacle_y, 
                             obstacle_r, low_speed_pos_x, low_speed_pos_y, 
                             low_speed_r);

    // store timing results
    for (int i = 0; i < tt[0].size(); i++){
        solve_timings_casadi_1.push_back(tt[0][i]);
        total_timings_casadi_1.push_back(tt[1][i]);
        update_timings_casadi_1.push_back(tt[1][i] - tt[0][i]);
    }

    // perform case CasADi Opti 2 (naive)
    tt = helper.performCasadiLoopReformulation(blocks, plotter, x0, 
                            corridor_1, corridor_2, nx, nu, nb_steps, 
                            nb_iterations, dt, T, nb_intervals, obstacle_x, 
                            obstacle_y, obstacle_r, low_speed_pos_x, 
                            low_speed_pos_y, low_speed_r);

    // store timing results
    for (int i = 0; i < tt[0].size(); i++){
        solve_timings_casadi_2.push_back(tt[0][i]);
        total_timings_casadi_2.push_back(tt[1][i]);
        update_timings_casadi_2.push_back(tt[1][i] - tt[0][i]);
    }

    // If only a single run is performed, write all the data to files
    if (nb_runs == 1){
        plotter.plotAnimationFramePython(nb_iterations-1);
        plotter.plotComputationTimesFullFull();
    }
    // else, write all timings results to files if this is the last run
    if (run_counter == nb_runs - 1){

        std::ofstream file_0("../../examples/plotting_data/timings_0.csv");
        if (file_0.is_open()){
            file_0<<"t_solve, t_update, t_total\n";
            for (int i = 0; i < solve_timings_casadi_1.size(); i++){
                file_0<<solve_timings_casadi_1[i]<<","<<
                    update_timings_casadi_1[i]<<","<<
                    total_timings_casadi_1[i]<<"\n";
            }
            file_0<<plotter.getMedian(solve_timings_casadi_1)<<","<<
                plotter.getMedian(update_timings_casadi_1)<<","<<
                plotter.getMedian(total_timings_casadi_1)<<"\n";

            file_0.close();
        }

        std::ofstream file_1("../../examples/plotting_data/timings_1.csv");
        if (file_1.is_open()){
            file_1<<"t_solve, t_update, t_total\n";
            for (int i = 0; i < solve_timings_casadi_2.size(); i++){
                file_1<<solve_timings_casadi_2[i]<<","<<
                    update_timings_casadi_2[i]<<","<<
                    total_timings_casadi_2[i]<<"\n";
            }
            file_1<<plotter.getMedian(solve_timings_casadi_2)<<","<<
                plotter.getMedian(update_timings_casadi_2)<<","<<
                plotter.getMedian(total_timings_casadi_2)<<"\n";

            file_1.close();
        } 

        std::ofstream file_2("../../examples/plotting_data/timings_2.csv");
        if (file_2.is_open()){
            file_2<<"t_solve, t_update, t_total\n";
            for (int i = 0; i < solve_timings_adaptive.size(); i++){
                file_2<<solve_timings_adaptive[i]<<","<<
                    update_timings_adaptive[i]<<","<<
                    total_timings_adaptive[i]<<"\n";
            }
            file_2<<plotter.getMedian(solve_timings_adaptive)<<","<<
                plotter.getMedian(update_timings_adaptive)<<","<<
                plotter.getMedian(total_timings_adaptive)<<"\n";

            file_2.close();
        }
    }
    }

    if (COMPUTE_SCALING_RESULTS){
        plotter.writeScaling();
    }

    }


    return 0;
}