#include "core.hpp"
#include "makeBuildingBlocks.hpp"
#include "plotter.hpp"
#include "adaptiveCorridorExampleHelpers.hpp"
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

int main(){
    //////////////////////////////////////
    // general definition of parameters //
    //////////////////////////////////////

    double T = 5.0;
    int Nmax = 100;
    int max_nb_extra_instances = 10;
    int nb_intervals = 8;
    int nb_reduced_horizon = 0;
    int nb_steps = 3;
    int nb_iterations = 130;
    std::vector<double> x0 = {-9.8, 1.0, 0.1, -0.1};
    double dt = T/((nb_intervals)*(nb_steps-1));
    double T_reduced = T - nb_reduced_horizon*(nb_steps-1)*dt;
    bool USE_PLOTTER = true;

    double view_radius = 10;
    std::vector<double> obstacle_x = 
        {40.0, 35.0, 45.0, 47.0, 46.0, 55.0, 50.5, 49};
    std::vector<double> obstacle_y = 
        {2.5,  1.2,  9.5,  7.0,  0.0,  8.5,  3.0,  2.5};
    std::vector<double> obstacle_r = 
        {2.0,  1.0,  3.3,  1.2,  1.9,  4.8,  0.8,  1.3};
    std::vector<double> low_speed_pos_x = 
        {3.9, 14, 23};
    std::vector<double> low_speed_pos_y = 
        {1.2, 3.0, 5.2};
    std::vector<double> low_speed_r = 
        {2.7, 2.7, 2.7};
    double low_speed_min = 0.0;

    ///////////////////////////////
    // Prepare for multiple runs //
    ///////////////////////////////
    int nb_runs = 1;
    std::vector<std::vector<double>> tt; // latest timing results
    std::vector<double> solve_timings_adaptive = {};
    std::vector<double> update_timings_adaptive = {};
    std::vector<double> total_timings_adaptive = {};
    std::vector<double> solve_timings_casadi_1 = {};
    std::vector<double> update_timings_casadi_1 = {};
    std::vector<double> total_timings_casadi_1 = {};
    std::vector<double> solve_timings_casadi_2 = {};
    std::vector<double> update_timings_casadi_2 = {};
    std::vector<double> total_timings_casadi_2 = {};

    for (int run_counter = 0; run_counter < nb_runs; run_counter++){
    cout<<endl<<endl<<endl<<endl;
    cout<<"STARTING RUN "<<run_counter<<endl;
    cout<<endl<<endl<<endl<<endl;

    /////////////////////////////
    // create object instances //
    /////////////////////////////

    std::vector<double> corridor_1;
    std::vector<double> corridor_2;
    BuildingBlocks blocks = makeBlocksCorridors(nb_steps, corridor_1, 
                                                corridor_2);
    int nx = blocks.get_nx();
    int nu = blocks.get_nu();

    // create plotter object
    Plotter plotter = Plotter(blocks.getFreeTime(), blocks.get_nx(), 
                              blocks.get_nu(), view_radius);
    std::string file_path = __FILE__;
    std::string dir_path = file_path.substr(0, file_path.rfind("/"));
    plotter.setOutputFolder(dir_path + "/plotting_data/");
    Plotter* plotter_ptr = nullptr;
    if (USE_PLOTTER){
        plotter_ptr = &plotter;
        for (int i = 0; i < obstacle_x.size(); i++){
            plotter.addCircle(obstacle_x[i], obstacle_y[i], obstacle_r[i]);
            plotter.setCircleVisibility(i, false);
        }
        for (int i = 0; i < low_speed_pos_x.size(); i++){
            plotter.addCircle(low_speed_pos_x[i], low_speed_pos_y[i], low_speed_r[i]);
            plotter.setCircleVisibility(obstacle_x.size() + i, false);
            plotter.setCircleColor(obstacle_x.size() + i, "b");
        }
        plotter.addCorridor(corridor_1);
        plotter.addCorridor(corridor_2);
    }
    
    // create adaptive NLP object
    // for (int nb_intervals : {8, 15, 25, 35, 45}){
    // double dt = T/((nb_intervals)*(nb_steps-1));
    // double T_reduced = T - nb_reduced_horizon*(nb_steps-1)*dt;
    // for (int Nmax : {30, 40, 50, 60, 75, 100, 120, 150, 200, 250, 300}){
    AdaptiveNLP adaptiveNLP = AdaptiveNLP(blocks, T_reduced, Nmax, 
                                          max_nb_extra_instances);
    std::vector<int> nks(nb_intervals-nb_reduced_horizon, nb_steps);
    adaptiveNLP.initTimeSteps(nks);
    std::vector<int> all(adaptiveNLP.getN());
    for (int i = 0; i < all.size(); i++){
        all[i] = i + (i >= adaptiveNLP.getFinalInd());
    }
    adaptiveNLP.addExtraConstraint(all, {corridor_1});

    // create helper object
    AdaptiveCorridorHelper helper = AdaptiveCorridorHelper(adaptiveNLP, 
                                                           view_radius, 
                                                           low_speed_min);

    // adaptive approach
    tt = helper.performAdaptiveLoop(USE_PLOTTER, adaptiveNLP, plotter, x0, 
                               corridor_1, corridor_2, nx, nu, nb_steps, 
                               nb_iterations, dt, nb_intervals, obstacle_x, 
                               obstacle_y, obstacle_r, low_speed_pos_x, 
                               low_speed_pos_y, low_speed_r);
    for (int i = 0; i < tt[0].size(); i++){
        solve_timings_adaptive.push_back(tt[0][i]);
        total_timings_adaptive.push_back(tt[1][i]);
        update_timings_adaptive.push_back(tt[1][i] - tt[0][i]);
    }

    // casadi::opti reference
    tt = helper.performCasadiLoop(USE_PLOTTER, blocks, plotter, x0, corridor_1, 
                             corridor_2, nx, nu, nb_steps, nb_iterations, dt, 
                             T, nb_intervals, obstacle_x, obstacle_y, 
                             obstacle_r, low_speed_pos_x, low_speed_pos_y, 
                             low_speed_r);
    for (int i = 0; i < tt[0].size(); i++){
        solve_timings_casadi_1.push_back(tt[0][i]);
        total_timings_casadi_1.push_back(tt[1][i]);
        update_timings_casadi_1.push_back(tt[1][i] - tt[0][i]);
    }

    // casadi::opti with reformulation (naive)
    tt = helper.performCasadiLoopReformulation(USE_PLOTTER, blocks, plotter, x0, 
                            corridor_1, corridor_2, nx, nu, nb_steps, 
                            nb_iterations, dt, T, nb_intervals, obstacle_x, 
                            obstacle_y, obstacle_r, low_speed_pos_x, 
                            low_speed_pos_y, low_speed_r);
    for (int i = 0; i < tt[0].size(); i++){
        solve_timings_casadi_2.push_back(tt[0][i]);
        total_timings_casadi_2.push_back(tt[1][i]);
        update_timings_casadi_2.push_back(tt[1][i] - tt[0][i]);
    }
    // }

    // // plotter.writeScaling();
    if (nb_runs == 1){
        plotter.plotAnimationFramePython(nb_iterations-1);
        plotter.plotComputationTimesFullFull();
    }

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

    // // plotter.plotIterCounts();
    }


    return 0;
}