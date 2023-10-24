#include "core.hpp"
#include "plotter.hpp"
#include "manyObstaclesHelper.hpp"
#include "makeBuildingBlocks.hpp"
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

int main(){
    //////////////////////////////////////
    // general definition of parameters //
    //////////////////////////////////////

    double T = 1.0;
    int Nmax = 60;
    int max_nb_extra_instances = 17;
    int nb_intervals = 30;
    int nb_reduced_horizon = 0;
    int nb_steps = 3;
    std::vector<double> start = {0.0, 0.0, 0.0, 0.0};
    std::vector<double> stop = {3.0, 0.0, 0.0};
    double dt = T/((nb_intervals)*(nb_steps-1));
    
    std::vector<double> obstacle_x_tot = 
        {1.5,  1.8,  0.5,  2.2,  1.4,  2.1,  2.5,  0.8,  2.9,  3.1, -0.2,  0.5,  1.5,  2.6,  1.9,  3.3, -0.4,  2.0,  1.0,  0.8,  2.8,  2.4};
    std::vector<double> obstacle_y_tot = 
        {-10, 20,  21,  22,  20, -20, -21, -22, -25,  30,  26,  35, -40, -50,  45,  48, -36,  29, -18,  50, -49,  27};
    std::vector<double> obstacle_r_tot = 
        {0.3,  0.4,  0.3,  0.2,  0.1,  0.1, 0.15,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.2 , 0.3,  0.4,  0.5,  0.6};
    
    /////////////////////////////
    // create object instances //
    /////////////////////////////

    BuildingBlocks blocks = makeBlocksManyObstacles();
    int nx = blocks.get_nx();
    int nu = blocks.get_nu();

    std::vector<double> timings_scaling_adaptive;
    std::vector<double> timings_scaling_casadi_1;
    std::vector<double> timings_scaling_casadi_2;
    std::vector<double> t_comp_scaling_casadi_2;
    std::vector<int> scaling_nb_obs;

    // create plotter object
    Plotter plotter = Plotter(blocks.getFreeTime(), blocks.get_nx(), 
                              blocks.get_nu());
    Plotter* plotter_ptr = nullptr;
    plotter_ptr = &plotter;
    for (int i = 0; i < obstacle_x_tot.size(); i++){
        plotter.addCircle(obstacle_x_tot[i], obstacle_y_tot[i], obstacle_r_tot[i]);
        plotter.setCircleVisibility(i, true);
    }

    for (int nb_runs = 0; nb_runs < 1; nb_runs++){
    for (int nb_obs = 1; nb_obs < obstacle_x_tot.size(); nb_obs++){

    scaling_nb_obs.push_back(nb_obs-1);
    std::vector<double> obstacle_x = std::vector<double>(nb_obs);
    std::vector<double> obstacle_y = std::vector<double>(nb_obs);
    std::vector<double> obstacle_r = std::vector<double>(nb_obs);
    for (int i = 0; i < nb_obs; i++){
        obstacle_x[i] = obstacle_x_tot[i];
        obstacle_y[i] = obstacle_y_tot[i];
        obstacle_r[i] = obstacle_r_tot[i];
    }

    // create adaptive NLP object
    AdaptiveNLP adaptiveNLP = AdaptiveNLP(blocks, T, Nmax, 
                                          max_nb_extra_instances);
    std::vector<int> nks(nb_intervals-nb_reduced_horizon, nb_steps);
    adaptiveNLP.initTimeSteps(nks);
    std::vector<int> all(adaptiveNLP.getN());
    for (int i = 0; i < all.size(); i++){
        all[i] = i + (i >= adaptiveNLP.getFinalInd());
    }

    // create helper object
    ManyObstaclesHelper helper = ManyObstaclesHelper(adaptiveNLP);

    // adaptive approach
    int max_nb_iterations = 10;
    std::vector<double> timings_adaptive(max_nb_iterations);
    std::vector<std::vector<std::vector<double>>> formatted_solutions_adaptive(
        max_nb_iterations, 
        std::vector<std::vector<double>>(nx+nu));
    std::vector<int> nb_constraints(max_nb_iterations);
    double total_time_adaptive;
    int nb_iterations_adaptive;
    // std::cout<<"Starting adaptive approach"<<std::endl;
    helper.performAdaptiveLoop(plotter, start, stop, nx, 
                               nu, obstacle_x, obstacle_y, 
                               obstacle_r, 
                               max_nb_iterations, 
                               timings_adaptive, 
                               formatted_solutions_adaptive, 
                               nb_constraints,
                               total_time_adaptive,
                               nb_iterations_adaptive);
    timings_scaling_adaptive.push_back(total_time_adaptive);

    // casadi::opti reference
    double total_time_casadi_1;
    std::vector<std::vector<std::vector<double>>> formatted_solutions_casadi_1(
        1, std::vector<std::vector<double>>(nx+nu));
    // std::cout<<"Starting casadi 1 approach"<<std::endl;
    helper.solveCompleteNLPCasadi(blocks, plotter, start, stop, nx, nu, 
                                  nb_steps, dt, nb_intervals, obstacle_x, 
                                  obstacle_y, obstacle_r, 
                                  formatted_solutions_casadi_1,
                                  total_time_casadi_1);
    timings_scaling_casadi_1.push_back(total_time_casadi_1);

    // casadi::opti with reformulation (naive)
    double total_time_casadi_2;
    int nb_iterations_casadi_2;
    std::vector<double> timings_casadi_2(max_nb_iterations);
    std::vector<std::vector<std::vector<double>>> formatted_solutions_casadi_2(
        max_nb_iterations, std::vector<std::vector<double>>(nx+nu));
    // std::cout<<"Starting casadi 2 approach"<<std::endl;
    helper.performCasadiLoopReformulation(blocks, plotter, start, stop, nx, nu,
                                          nb_steps, dt, nb_intervals, 
                                          obstacle_x, obstacle_y, obstacle_r,
                                          max_nb_iterations,
                                          timings_casadi_2,
                                          formatted_solutions_casadi_2,
                                          total_time_casadi_2,
                                          nb_iterations_casadi_2);
    timings_scaling_casadi_2.push_back(total_time_casadi_2);
    t_comp_scaling_casadi_2.push_back(sum(timings_casadi_2));

    // std::cout<<"Total time adaptive approach: "<<total_time_adaptive<<" ms (";
    // for (int i = 0; i < nb_iterations_adaptive; i++){
    //     cout<<timings_adaptive[i]<<" ";
    // }
    // cout<<")"<<endl;
    // std::cout<<"Total time casadi 1 appraoch: "<<total_time_casadi_1<<" ms"<<std::endl;
    // std::cout<<"Total time casadi 2 approach: "<<total_time_casadi_2<<" ms (";
    // for (int i = 0; i < nb_iterations_casadi_2; i++){
    //     cout<<timings_casadi_2[i]<<" ";
    // }
    // cout<<")"<<endl;

    // if (nb_obs == obstacle_x_tot.size()-1){
    //     for (int i = 0; i < nb_iterations_adaptive; i++){
    //         plotter.plotTrajectory(formatted_solutions_adaptive[i][0], 
    //                             formatted_solutions_adaptive[i][1]);
    //     }
    //     plotter.plotTrajectory(formatted_solutions_casadi_1[0][0],
    //                         formatted_solutions_casadi_1[0][1]);
    // 
    // }
    
    }
    }

    // std::cout<<"adaptive: "; for (auto e : timings_scaling_adaptive){std::cout<<e<<" ";}std::cout<<endl;
    // std::cout<<"casadi 1: "; for (auto e : timings_scaling_casadi_1){std::cout<<e<<" ";}std::cout<<endl;
    // std::cout<<"casadi 2: "; for (auto e : timings_scaling_casadi_2){std::cout<<e<<" ";}std::cout<<endl;

    std::ofstream file("../../examples/plotting_data/virtualObstaclesScaling.csv");
    if (file.is_open()){
        file << "nb_obs, adaptive, casadi 1, casadi 2\n";
        for (int i = 0; i < timings_scaling_adaptive.size(); i++){
            file << scaling_nb_obs[i]<<","<<timings_scaling_adaptive[i]<<","
                <<timings_scaling_casadi_1[i]<<
                ","<<timings_scaling_casadi_2[i]<<", "<<
                t_comp_scaling_casadi_2[i]<<"\n";
        }
        file.close();
    }

    return 0;
}