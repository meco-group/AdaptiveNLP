#include "../include/plotter.hpp"
#include <vector>
#include <optional>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <cmath>

// #include "matplotlibcpp.h"
// namespace plt = matplotlibcpp;

Plotter::Plotter(bool free_time, int nx, int nu, double view_radius){
    free_time_ = free_time;
    nx_ = nx;
    nu_ = nu;
    view_radius_ = view_radius;

    xk_ = std::vector<double>(nx);
    uk_ = std::vector<double>(nu);

    scaling_median_total_times_adaptive_ = {};
    scaling_median_total_times_ref_ = {};
    scaling_median_total_times_naive_ = {};
    
    scaling_median_times_adaptive_ = {};
    scaling_median_times_ref_ = {};
    scaling_median_times_naive_ = {};
    
    scaling_Ns_adaptive_ = {};
    scaling_Ns_ref_ = {};
    scaling_Ns_naive_ = {};
}

Plotter::Plotter(bool free_time, int nx, int nu, int N, int final_ind,
                std::vector<std::optional<int>>& next_ind,
                std::vector<std::optional<double>>& time_from_ind){
    free_time_ = free_time;
    nx_ = nx;
    nu_ = nu;
    N_ = N;
    final_ind_ = final_ind;
    next_ind_ = std::vector<std::optional<int>>(N_);
    time_x_ = std::vector<double>(N_+1);
    time_u_ = std::vector<double>(N_);
    int i = 0;
    int time_ptr = 0;
    while (i != final_ind){
        next_ind_[i] = next_ind[i].value();
        time_x_[time_ptr] = time_from_ind[i].value();
        time_u_[time_ptr] = time_from_ind[i].value();
        time_ptr++;
        i = next_ind[i].value();
    }
    xk_ = std::vector<double>(nx);
    uk_ = std::vector<double>(nu);
}

void Plotter::update(int N, int final_ind,
                     std::vector<std::optional<int>>& next_ind,
                     std::vector<std::optional<double>>& time_from_ind){
    N_ = N;
    final_ind_ = final_ind;
    next_ind_ = std::vector<std::optional<int>>(N_+1);
    time_x_ = std::vector<double>(N_+1);
    time_u_ = std::vector<double>(N_);
    int i = 0;
    int time_ptr = 0;
    while (i != final_ind){
        next_ind_[i] = next_ind[i].value();
        time_x_[time_ptr] = time_from_ind[i].value();
        time_u_[time_ptr] = time_from_ind[i].value();
        time_ptr++;
        i = next_ind[i].value();
    }
};

// void Plotter::plotTrajectory(const double* sol, std::string c, 
//                              bool new_figure){
//     std::vector<double> xx(N_+1);
//     std::vector<double> yy(N_+1);

//     int k = 0;
//     for (int i = 0; i < N_+1; i++){
//         get_x(sol, k);
//         xx[i] = xk_[0];
//         yy[i] = xk_[1];
//         if (k != final_ind_){
//             k = next_ind_[k].value();
//         }

//     }
//     return plotTrajectory(xx, yy, c, new_figure);
// }

// void Plotter::plotTrajectory(std::vector<double>& xx, std::vector<double>& yy,
//                              std::string c, bool new_figure){
//     if (new_figure){ plt::figure();}
//     for (int i = 0; i < corridors_.size(); i++){
//         plotCorridor(i);
//     }
//     for (int i = 0; i < circle_x_.size(); i++){
//         plotCircle(i);
//     }
//     plt::plot(xx, yy, {{"color", c}, {"marker", "o"}});
//     plt::axis("equal");
//     if (new_figure){
//         plt::save("figures/trajectory_" + std::to_string(plotter_count_) + 
//             file_format_, dpi_);
//         plt::close();
//     }
//     plotter_count_++;
// }

// void Plotter::plotControls(const double* sol){
//     std::vector<std::vector<double>> uu(nu_, std::vector<double>(N_, 0.0));
//     int k = 0;
//     for (int i = 0; i < N_; i++){
//         get_u(sol, k);
//         for (int j = 0; j < nu_; j++){
//             uu[j][i] = uk_[j];
//         }
//         k = next_ind_[k].value();
//     }
//     return plotControls(uu);
// };

// void Plotter::plotControls(std::vector<std::vector<double>>& uu){
//     std::vector<std::string> labels = {"a", "delta"};
//     plt::figure();
//     for (int i = 0; i < nu_; i++){
//         plt::plot(time_u_, uu[i], {{"marker", "o"}, {"label", labels[i]}});
//     }
//     plt::legend();
//     // plt.xlabel()
//     plt::save("figures/controls_" + std::to_string(plotter_count_) + 
//               file_format_, dpi_);
//     plt::close();
//     plotter_count_++;
// }

void Plotter::writeControlsToFile(const double* sol){
    std::vector<std::vector<double>> uu(nu_, std::vector<double>(N_, 0.0));
    int k = 0;
    for (int i = 0; i < N_; i++){
        get_u(sol, k);
        for (int j = 0; j < nu_; j++){
            uu[j][i] = uk_[j];
        }
        k = next_ind_[k].value();
    }
    return writeControlsToFile(uu);
};

void Plotter::writeControlsToFile(std::vector<std::vector<double>>& uu){
    std::ofstream file("plotting_data/controls.csv");
    if (file.is_open()){
        file<<"t";
        for (int i = 0; i < nu_; i++){
            file<<", u_"<<i;
        }
        file<<"\n";

        for (int i = 0; i < N_; i++){
            file<<time_u_[i];
            for (int j = 0; j < nu_; j++){
                file<<", "<<uu[j][i];
            }
            file<<"\n";
        }
        file.close();
    }
}; 

// void Plotter::plotTravelledTrajectory(std::vector<std::vector<double>> x0s, 
//                                       std::string c, int frame_nb){
//     std::vector<double> xx(frame_nb+1);
//     std::vector<double> yy(frame_nb+1);
//     for (int i = 0; i < frame_nb+1; i++){
//         xx[i] = x0s[i][0]; yy[i] = x0s[i][1];
//     }
//     plt::plot(xx, yy, {{"color", c}});
// };

// void Plotter::plotAnimationFrame(int frame_nb){
//     if (!ref_solution_available_ && !adaptive_solution_available_ &&
//             !naive_solution_available_){
//         return;
//     }
//     // plot trajectories
//     plt::subplot2grid(3, 1, 0, 0);
//     plt::cla();
//     if (adaptive_solution_available_){
//         for (int i = 0; i < visible_circles_.size(); i++){
//             visible_circles_[i] = visible_constraints_adaptive_[frame_nb][i];
//         }
//     } else {
//         for (int i = 0; i < visible_circles_.size(); i++){
//             visible_circles_[i] = visible_constraints_ref_[frame_nb][i];
//         }
//     }
//     if (naive_solution_available_){
//         plotTravelledTrajectory(x0s_naive_, "orange", frame_nb);
//         plotTrajectory(formatted_solutions_naive_[frame_nb][0], 
//                        formatted_solutions_naive_[frame_nb][1], "orange", false);
//         plotCircle(x0s_naive_[frame_nb][0], x0s_naive_[frame_nb][1], 
//                    view_radius_, "k");
//     }
//     if (ref_solution_available_){
//         plotTravelledTrajectory(x0s_ref_, "0.3", frame_nb);
//         plotTrajectory(formatted_solutions_ref_[frame_nb][0], 
//                        formatted_solutions_ref_[frame_nb][1], "0.3", false);
//         plotCircle(x0s_ref_[frame_nb][0], x0s_ref_[frame_nb][1], 
//                    view_radius_, "k");
//     }
//     if (adaptive_solution_available_){
//         plotTravelledTrajectory(x0s_adaptive_, "k", frame_nb);
//         plotTrajectory(formatted_solutions_adaptive_[frame_nb][0], 
//                        formatted_solutions_adaptive_[frame_nb][1], "b", false);
//         plotCircle(x0s_adaptive_[frame_nb][0], x0s_adaptive_[frame_nb][1], 
//                    view_radius_, "k");
//     }
        
//     double veh_width = 0.8; double veh_height = 1.5;
//     double tilt;
//     std::vector<double> center;
//     if (adaptive_solution_available_){
//         tilt = x0s_adaptive_[frame_nb][3];
//         center = {x0s_adaptive_[frame_nb][0], 
//                   x0s_adaptive_[frame_nb][1]};
//     } else {
//         tilt = x0s_ref_[frame_nb][3];
//         center = {x0s_ref_[frame_nb][0], 
//                   x0s_ref_[frame_nb][1]};
//     }

//     std::vector<double> veh_x = {center[0] - veh_height/2*cos(tilt) + 
//                                              veh_width/2*sin(tilt),
//                                  center[0] + veh_height/2*cos(tilt) + 
//                                              veh_width/2*sin(tilt),
//                                  center[0] + veh_height/2*cos(tilt) - 
//                                              veh_width/2*sin(tilt),
//                                  center[0] - veh_height/2*cos(tilt) - 
//                                              veh_width/2*sin(tilt)};
//     std::vector<double> veh_y = {center[1] - veh_width/2*cos(tilt) - 
//                                              veh_height/2*sin(tilt),
//                                  center[1] - veh_width/2*cos(tilt) + 
//                                              veh_height/2*sin(tilt),
//                                  center[1] + veh_width/2*cos(tilt) + 
//                                              veh_height/2*sin(tilt),
//                                  center[1] + veh_width/2*cos(tilt) - 
//                                              veh_height/2*sin(tilt)};
//     plotVehicle(veh_x,veh_y);
//     plt::axis("equal");
//     plt::xlim(-15, 65);
//     plt::ylim(-5,12);

//     // plot number of constraints
//     plt::subplot2grid(3, 1, 1, 0);
//     plt::cla();
//     if (naive_solution_available_){
//         std::vector<int> nb_naive(frame_nb+1);
//         for (int i = 0; i < frame_nb+1; i++){
//             nb_naive[i] = nb_constraints_naive_[i];
//         }
//         plt::plot(nb_naive, {{"color", "orange"}, {"linewidth", "3.5"}});
//         plt::xlim(0, int(nb_constraints_naive_.size()));
//         plt::ylim(
//             *std::min_element(nb_constraints_naive_.begin(), 
//                               nb_constraints_naive_.end()) - 10,
//             *std::max_element(nb_constraints_naive_.begin(), 
//                               nb_constraints_naive_.end()) + 10);
//     }
//     if (ref_solution_available_){
//         std::vector<int> nb_ref(frame_nb+1);
//         for (int i = 0; i < frame_nb+1; i++){
//             nb_ref[i] = nb_constraints_ref_[i];
//         }
//         plt::plot(nb_ref, {{"color", "0.3"}});
//         plt::xlim(0, int(nb_constraints_ref_.size()));
//         plt::ylim(
//             *std::min_element(nb_constraints_ref_.begin(), 
//                               nb_constraints_ref_.end()) - 10,
//             *std::max_element(nb_constraints_ref_.begin(), 
//                               nb_constraints_ref_.end()) + 10);
//     }
//     if (adaptive_solution_available_){
//         std::vector<int> nb_adaptive(frame_nb+1);
//         std::vector<std::vector<int>> nb_adaptive_specific(
//             nb_constraints_specific_adaptive_[0].size() - 2, 
//             std::vector<int>(frame_nb+1));
//         for (int i = 0; i < frame_nb+1; i++){
//             nb_adaptive[i] = nb_constraints_adaptive_[i];
//             for (int j = 2; j < nb_constraints_specific_adaptive_[0].size(); j++){
//                 nb_adaptive_specific[j-2][i] = 
//                     nb_constraints_specific_adaptive_[i][j];
//             }
//         }
//         assert (nb_adaptive_specific.size() < 5);
//         std::vector<std::string> colors = {"blue", "dodgerblue", "navy", 
//                                            "mediumslateblue", "purple"};
//         // for (int j = 2; j < nb_constraints_specific_adaptive_[0].size(); j++){
//         for (int j = nb_constraints_specific_adaptive_[0].size()-1; j >= 2; 
//                 j--){
//             plt::plot(nb_adaptive_specific[j-2], {{"color", colors[j-2]}, 
//                                                 // {"linestyle", ":"},
//                                                 // {"linewidth", "1.0"}
//                                                 });
//         }
//         // plt::plot(nb_adaptive, {{"color", "b"}, {"linewidth", "1.3"}});
//         plt::xlim(0.0, double(nb_constraints_adaptive_.size()));
//         if (!ref_solution_available_){
//             plt::ylim(//0,
//                 *std::min_element(nb_constraints_adaptive_.begin(), 
//                                   nb_constraints_adaptive_.end()) - 10,
//                 *std::max_element(nb_constraints_adaptive_.begin(), 
//                                   nb_constraints_adaptive_.end()) + 10);
//         } else {
//             plt::ylim(
//                     // 0,
//                       std::min({
//                         *std::min_element(nb_constraints_adaptive_.begin(), 
//                                           nb_constraints_adaptive_.end()),
//                         *std::min_element(nb_constraints_ref_.begin(), 
//                                           nb_constraints_ref_.end())}) - 10, 
//                       std::max({
//                         *std::max_element(nb_constraints_adaptive_.begin(), 
//                                           nb_constraints_adaptive_.end()),
//                         *std::max_element(nb_constraints_ref_.begin(), 
//                                           nb_constraints_ref_.end())}) + 10);
//         }
//     }
//     plt::ylabel("number of\nconstraints");

//     // plot computation times
//     plt::subplot2grid(3, 1, 2, 0);
//     plt::cla();
//     if (naive_solution_available_){
//         std::vector<int> times_naive(frame_nb+1);
//         std::vector<int> total_times_naive(frame_nb+1);
//         for (int i = 0; i < frame_nb+1; i++){
//             total_times_naive[i] = total_timings_naive_[i];
//             times_naive[i] = timings_naive_[i];
//         }
//         plt::plot(times_naive, {{"color", "orange"}, {"linewidth", "0.4"}});
//         plt::plot(total_times_naive, {{"color", "orange"}});
//         plt::xlim(0, int(timings_naive_.size()));
//         plt::ylim(0.0, *std::max_element(total_timings_naive_.begin(), 
//                                          total_timings_naive_.end()) + 5);
//     }
//     if (ref_solution_available_){
//         std::vector<int> times_ref(frame_nb+1);
//         std::vector<int> total_times_ref(frame_nb+1);
//         for (int i = 0; i < frame_nb+1; i++){
//             total_times_ref[i] = total_timings_ref_[i];
//             times_ref[i] = timings_ref_[i];
//         }
//         plt::plot(times_ref, {{"color", "0.3"}, {"linewidth", "0.4"}});
//         plt::plot(total_times_ref, {{"color", "0.3"}});
//         plt::xlim(0, int(timings_ref_.size()));
//         plt::ylim(0.0, *std::max_element(total_timings_ref_.begin(), 
//                                          total_timings_ref_.end()) + 5);
//     }
//     if (adaptive_solution_available_){
//         std::vector<int> times_adaptive(frame_nb+1);
//         std::vector<int> total_times_adaptive(frame_nb+1);
//         for (int i = 0; i < frame_nb+1; i++){
//             total_times_adaptive[i] = total_timings_adaptive_[i];
//             times_adaptive[i] = timings_adaptive_[i];
//         }
//         plt::plot(times_adaptive, {{"color", "b"}, {"linewidth", "0.4"}});
//         plt::plot(total_times_adaptive, {{"color", "b"}});
//         plt::xlim(0, int(timings_adaptive_.size()));
//         if (!ref_solution_available_){
//             plt::ylim(0.0, *std::max_element(total_timings_adaptive_.begin(), 
//                                              total_timings_adaptive_.end()) + 5);
//         } else if (naive_solution_available_){
//             plt::ylim(0.0, std::max({
//                         *std::max_element(total_timings_adaptive_.begin(), 
//                                           total_timings_adaptive_.end()),
//                         *std::max_element(total_timings_ref_.begin(), 
//                                           total_timings_ref_.end()),
//                         *std::max_element(total_timings_naive_.begin(), 
//                                           total_timings_naive_.end())}) + 5);
//         } else {
//             plt::ylim(0.0, std::max({
//                         *std::max_element(total_timings_adaptive_.begin(), 
//                                           total_timings_adaptive_.end()),
//                         *std::max_element(total_timings_ref_.begin(), 
//                                           total_timings_ref_.end())}) + 5);
//         }
//     }

//     plt::ylabel("computation\ntime [ms]");
//     plt::xlabel("iteration number");

//     plt::save("figures/animation/animation_frame_"+std::to_string(frame_nb)+
//         file_format_, dpi_);
// };

template <typename T>
void writeRowToFile(std::ofstream& file, std::vector<T>& row){
    for (int i = 0; i < row.size(); i++){
        file << row[i];
        if (i < row.size()-1){
            file << ",";
        }
    }
    file <<"\n";
}

void Plotter::plotAnimationFramePython(int frame_nb){
    if (!ref_solution_available_ && !adaptive_solution_available_ &&
            !naive_solution_available_){
        return;
    }

    // write circle information
    std::ofstream file1("plotting_data/circles.csv");
    if (file1.is_open()){
        for (int i = 0; i < circle_x_.size(); i++){
            file1 <<"circle "<<i;
            if (i < circle_x_.size()-1){
                file1 <<",";
            }
        }
        file1<<"\n";
        writeRowToFile(file1, circle_x_);
        writeRowToFile(file1, circle_y_);
        writeRowToFile(file1, circle_r_);
        writeRowToFile(file1, circle_colors_);

        for (int i = 0; i < visible_constraints_adaptive_.size(); i++){
            writeRowToFile(file1, visible_constraints_adaptive_[i]);
        }

        file1.close();
    }

    // write corridor information
    std::ofstream file2("plotting_data/corridors.csv");
    if (file2.is_open()){
        file2 << "x_min, x_max, y_min, y_max\n"; 
        for (int i = 0; i < corridors_.size(); i++){
            writeRowToFile(file2, corridors_[i]);
        }
        file2.close();
    }

    // write simulation results
    if (adaptive_solution_available_){
        std::ofstream file3("plotting_data/results_adaptive.csv");
        if (file3.is_open()){
            file3 << "timings, timings_total, nb_constraints, nb_specific_1, "<<
                    "nb_specific_2, nb_specific_3, nb_specific_4, x0s[0], "<<
                    "x0s[1], x0s[2], x0s[3]\n";
            for (int i = 0; i < timings_adaptive_.size(); i++){
                file3 << timings_adaptive_[i]<<","
                    <<total_timings_adaptive_[i]<<
                    ","<<nb_constraints_adaptive_[i]<<","<<
                    nb_constraints_specific_adaptive_[i][2]<<","<<
                    nb_constraints_specific_adaptive_[i][3]<<","<<
                    nb_constraints_specific_adaptive_[i][4]<<","<<
                    nb_constraints_specific_adaptive_[i][5]<<","<<
                    x0s_adaptive_[i][0]<<","<<x0s_adaptive_[i][1]<<","<<
                    x0s_adaptive_[i][2]<<","<<x0s_adaptive_[i][3]<<"\n";
            }
            file3.close();
        }

        std::ofstream file4("plotting_data/solutions_adaptive.csv");
        for (int i = 0; i < formatted_solutions_adaptive_.size(); i++){
            for (int j = 0; j < formatted_solutions_adaptive_[i].size(); j++){
                file4<<formatted_solutions_adaptive_[i][j].size()<<"\n";
                for (int k = 0; k < formatted_solutions_adaptive_[i][j].size();
                        k++){
                    file4<<formatted_solutions_adaptive_[i][j][k]<<"\n";
                }
            }
        }
        file4.close();
    }
    if (ref_solution_available_){
        std::ofstream file5("plotting_data/results_ref.csv");
        if (file5.is_open()){
            file5 << "timings, timings_total, nb_constraints, x0s[0], "<<
                    "x0s[1], x0s[2], x0s[3]\n";
            for (int i = 0; i < timings_ref_.size(); i++){
                file5 << timings_ref_[i]<<","<<total_timings_ref_[i]<<
                    ","<<nb_constraints_ref_[i]<<","<<
                    x0s_ref_[i][0]<<","<<x0s_ref_[i][1]<<","<<
                    x0s_ref_[i][2]<<","<<x0s_ref_[i][3]<<"\n";
            }
            file5.close();
        }

        std::ofstream file6("plotting_data/solutions_ref.csv");
        for (int i = 0; i < formatted_solutions_ref_.size(); i++){
            for (int j = 0; j < formatted_solutions_ref_[i].size(); j++){
                file6<<formatted_solutions_ref_[i][j].size()<<"\n";
                for (int k = 0; k < formatted_solutions_ref_[i][j].size();
                        k++){
                    file6<<formatted_solutions_ref_[i][j][k]<<"\n";
                }
            }
        }
        file6.close();
    }
    if (naive_solution_available_){
        std::ofstream file7("plotting_data/results_naive.csv");
        if (file7.is_open()){
            file7 << "timings, timings_total, nb_constraints, x0s[0], "<<
                    "x0s[1], x0s[2], x0s[3]\n";
            for (int i = 0; i < timings_naive_.size(); i++){
                file7 << timings_naive_[i]<<","<<total_timings_naive_[i]<<
                    ","<<nb_constraints_naive_[i]<<","<<
                    x0s_naive_[i][0]<<","<<x0s_naive_[i][1]<<","<<
                    x0s_naive_[i][2]<<","<<x0s_naive_[i][3]<<"\n";
            }
            file7.close();
        }
        std::ofstream file8("plotting_data/solutions_naive.csv");
        for (int i = 0; i < formatted_solutions_naive_.size(); i++){
            for (int j = 0; j < formatted_solutions_naive_[i].size(); j++){
                file8<<formatted_solutions_naive_[i][j].size()<<"\n";
                for (int k = 0; k < formatted_solutions_naive_[i][j].size();
                        k++){
                    file8<<formatted_solutions_naive_[i][j][k]<<"\n";
                }
            }
        }
    }

    
};

// void Plotter::plotControlAnimationFrame(int frame_nb){
//     if (!adaptive_solution_available_ && !ref_solution_available_){
//         return;
//     }
//     // plt::figure();
//     std::vector<double> y_lower = {-1.0, -1.047};
//     std::vector<double> y_upper = { 1.0,  1.047};

//     for (int i = 0; i < nu_; i++){
//         plt::subplot2grid(nu_, 1, i, 0);
//         plt::cla();
//         plt::ylim(y_lower[i], y_upper[i]);
        
//         if (adaptive_solution_available_){
//             plt::plot(formatted_solutions_adaptive_[frame_nb][nx_+i],
//                       {{"color", "b"}});
//         }
//         if (ref_solution_available_){
//             plt::plot(formatted_solutions_ref_[frame_nb][nx_+i],
//                       {{"color", "0.3"}});
//             plt::plot(formatted_solutions_ref_[frame_nb][nx_+i]);
//         }
//     }
//     plt::save("figures/animation_controls/animation_frame_"+
//               std::to_string(frame_nb)+file_format_, dpi_);
//     plt::close();
// };

void Plotter::addAdaptiveResults(
            std::vector<std::vector<std::vector<double>>> 
                formatted_solutions_adaptive,
            std::vector<std::vector<double>> x0s_adaptive,
            std::vector<std::vector<bool>> visible_constraints_adaptive,
            std::vector<int> nb_constraints_adaptive,
            std::vector<std::vector<int>> nb_constraints_specific_adaptive,
            std::vector<double> timings_adaptive,
            std::vector<double> total_timings_adaptive){
    formatted_solutions_adaptive_ = formatted_solutions_adaptive;
    x0s_adaptive_ = x0s_adaptive;
    visible_constraints_adaptive_ = visible_constraints_adaptive;
    nb_constraints_adaptive_ = nb_constraints_adaptive;
    nb_constraints_specific_adaptive_ = nb_constraints_specific_adaptive;
    timings_adaptive_ = timings_adaptive;
    total_timings_adaptive_ = total_timings_adaptive;
    adaptive_solution_available_ = true;
    scaling_median_total_times_adaptive_.push_back(
        getMedian(total_timings_adaptive_));
    scaling_median_times_adaptive_.push_back(
        getMedian(timings_adaptive_));
    scaling_Ns_adaptive_.push_back(
        formatted_solutions_adaptive_[timings_adaptive.size()-1][0].size());
};

void Plotter::addAdaptiveIterCounts(std::vector<int>& iter_counts){
    iter_counts_adaptive_ = iter_counts;
    adaptive_iter_counts_available_ = true;
}

void Plotter::addRefResults(
            std::vector<std::vector<std::vector<double>>>
                formatted_solutions_ref,
            std::vector<std::vector<double>> x0s_ref,
            std::vector<std::vector<bool>> visible_constraints_ref,
            std::vector<int> nb_constraints_ref,
            std::vector<std::vector<int>> nb_constraints_specific_ref,
            std::vector<double> timings_ref,
            std::vector<double> total_timings_ref){
    formatted_solutions_ref_ = formatted_solutions_ref;
    x0s_ref_ = x0s_ref;
    visible_constraints_ref_ = visible_constraints_ref;
    nb_constraints_ref_ = nb_constraints_ref;
    nb_constraints_specific_ref_ = nb_constraints_specific_ref;
    timings_ref_ = timings_ref;
    total_timings_ref_ = total_timings_ref;
    ref_solution_available_ = true;
    scaling_median_total_times_ref_.push_back(getMedian(total_timings_ref_));
    scaling_median_times_ref_.push_back(getMedian(timings_ref_));
    scaling_Ns_ref_.push_back(
        formatted_solutions_ref_[timings_ref.size()-1][0].size());
};

void Plotter::addRefIterCounts(std::vector<int>& iter_counts){
    iter_counts_ref_ = iter_counts;
    ref_iter_counts_available_ = true;
}

void Plotter::addNaiveResults(
            std::vector<std::vector<std::vector<double>>>
                formatted_solutions_naive,
            std::vector<std::vector<double>> x0s_naive,
            std::vector<std::vector<bool>> visible_constraints_naive,
            std::vector<int> nb_constraints_naive,
            std::vector<std::vector<int>> nb_constraints_specific_naive,
            std::vector<double> timings_naive,
            std::vector<double> total_timings_naive){
    formatted_solutions_naive_ = formatted_solutions_naive;
    x0s_naive_ = x0s_naive;
    visible_constraints_naive_ = visible_constraints_naive;
    nb_constraints_naive_ = nb_constraints_naive;
    nb_constraints_specific_naive_ = nb_constraints_specific_naive;
    timings_naive_ = timings_naive;
    total_timings_naive_ = total_timings_naive;
    naive_solution_available_ = true;
    scaling_median_total_times_naive_.push_back(
        getMedian(total_timings_naive_));
    scaling_median_times_naive_.push_back(getMedian(timings_naive_));
    scaling_Ns_naive_.push_back(
        formatted_solutions_naive_[timings_naive.size()-1][0].size());
};

void Plotter::addNaiveIterCounts(std::vector<int>& iter_counts){
    iter_counts_naive_ = iter_counts;
    naive_iter_counts_available_ = true;
}

// void Plotter::plotComputationTimes(std::vector<double>& timings){
//     plt::figure();
//     // plt::boxplot(timings, {{"widths", {"0.6"}}});
//     plt::boxplot(timings);
//     std::vector<double> ticks = {1.0};
//     std::vector<std::string> labels = {"adaptiveNLP"};
//     plt::xticks(ticks, labels);
//     plt::ylabel("computation times [ms]");
//     plt::grid(true);
//     plt::save("figures/computation_times" + file_format_, dpi_);
//     plt::close();
// };

// void Plotter::plotComputationTimes(std::vector<std::vector<double>>& timings){
//     plt::figure();
//     plt::boxplot(timings);
//     plt::ylabel("computation times [ms]");
//     plt::save("figures/computation_times" + file_format_, dpi_);
//     plt::close();
// };

// void Plotter::plotComputationTimes(){
//     if (ref_solution_available_ && adaptive_solution_available_){
//         plt::figure();
//         std::vector<std::vector<double>> timings = {timings_adaptive_, 
//                                                     timings_ref_};
//         plt::boxplot(timings);
//         std::vector<double> ticks = {1.0, 2.0};
//         std::vector<std::string> labels = {"adaptiveNLP", "casadi::Opti"};
//         plt::xticks(ticks, labels);
//         plt::ylabel("computation times [ms]");
//         plt::grid(true);

//         double median_adaptive = getMedian(timings_adaptive_);
//         std::ostringstream oss_adaptive;
//         oss_adaptive<<std::fixed<<std::setprecision(2)<<std::setw(5)<<
//             std::round(median_adaptive*100.0)/100.0;
//         plt::text(ticks[0] + 0.1, median_adaptive, oss_adaptive.str());

//         double median_ref = getMedian(timings_ref_);
//         std::ostringstream oss_ref;
//         oss_ref<<std::fixed<<std::setprecision(2)<<std::setw(5)<<
//             std::round(median_ref*100.0)/100.0;
//         plt::text(ticks[1] + 0.1, median_ref, oss_ref.str());
//         plt::save("figures/computation_times" + file_format_, dpi_);
//         plt::close();
//     } else if (adaptive_solution_available_){
//         plotComputationTimes(timings_adaptive_);
//     }
// };

// void Plotter::plotComputationTimesFull(){
//     if (ref_solution_available_ && adaptive_solution_available_){
//         plt::figure();
//         std::vector<std::vector<double>> timings = {timings_adaptive_,
//                                                     total_timings_adaptive_,
//                                                     timings_ref_,
//                                                     total_timings_ref_};
//         // std::vector<std::string> positions = {"1.0", "1.3", "2.0", "2.3"};
//         // plt::boxplot(timings, positions);
//         plt::boxplot(timings);
//         // std::vector<double> ticks = {1.0, 1.3, 2.0, 2.3};
//         std::vector<double> ticks = {1.0, 2.0, 3.0, 4.0};
//         std::vector<std::string> labels = {"adaptiveNLP\nt_solve", 
//                                            "adaptiveNLP\nt_total",
//                                            "casadi::Opti\nt_solve",
//                                            "casadi::Opti\nt_total"};
//         plt::xticks(ticks, labels);
//         plt::xlim(0.3, 4.7);
//         plt::ylabel("computation times [ms]");
//         plt::grid(false);

//         double text_shift = 0.25;

//         double median_adaptive = getMedian(timings_adaptive_);
//         plt::text(ticks[0] + text_shift, median_adaptive,
//                   getMedianText(timings_adaptive_));

//         double median_adaptive_total = getMedian(total_timings_adaptive_);
//         plt::text(ticks[1] + text_shift, median_adaptive_total,
//                   getMedianText(total_timings_adaptive_));

//         double median_ref = getMedian(timings_ref_);
//         plt::text(ticks[2] + text_shift, median_ref,
//                   getMedianText(timings_ref_));

//         double median_ref_total = getMedian(total_timings_ref_);
//         plt::text(ticks[3] + text_shift, median_ref_total,
//                   getMedianText(total_timings_ref_));
        
//         plt::save("figures/computation_times_full" + file_format_, dpi_);
//         plt::close();
//     }
// };

void Plotter::plotComputationTimesFullFull(){
    double y_max = 0.0;
    std::vector<std::vector<double>> timings;
    std::vector<std::string> labels;
    if (ref_solution_available_){
        std::vector<double> diff_ref(timings_ref_.size());
        for (int i = 0; i < timings_ref_.size(); i++){
            diff_ref[i] = total_timings_ref_[i] - 
                timings_ref_[i];
        }
        timings.push_back(timings_ref_);
        timings.push_back(diff_ref);
        timings.push_back(total_timings_ref_);
        labels.push_back("$Casadi::Opti\nt_{solve}$");
        labels.push_back("$Casadi::Opti\nt_{update}$");
        labels.push_back("$Casadi::Opti\nt_{total}$");
        y_max = std::max({y_max, *std::max_element(total_timings_ref_.begin(), 
                                                total_timings_ref_.end()),});
    }
    if (naive_solution_available_){
        std::vector<double> diff_naive(timings_naive_.size());
        for (int i = 0; i < timings_naive_.size(); i++){
            diff_naive[i] = total_timings_naive_[i] - 
                timings_naive_[i];
        }
        timings.push_back(timings_naive_);
        timings.push_back(diff_naive);
        timings.push_back(total_timings_naive_);
        labels.push_back("$Casadi::Opti\nt_{solve}$");
        labels.push_back("$Casadi::Opti\nt_{update}$");
        labels.push_back("$Casadi::Opti\nt_{total}$");
        y_max = std::max({y_max, *std::max_element(total_timings_ref_.begin(), 
                                                total_timings_ref_.end()),});
    }
    if (adaptive_solution_available_){
        std::vector<double> diff_adaptive(timings_adaptive_.size());
        for (int i = 0; i < timings_adaptive_.size(); i++){
            diff_adaptive[i] = total_timings_adaptive_[i] - 
                timings_adaptive_[i];
        }
        timings.push_back(timings_adaptive_);
        timings.push_back(diff_adaptive);
        timings.push_back(total_timings_adaptive_);
        labels.push_back("$AdaptiveNLP\nt_{solve}$");
        labels.push_back("$AdaptiveNLP\nt_{update}$");
        labels.push_back("$AdaptiveNLP\nt_{total}$");
        y_max = std::max({y_max, *std::max_element(total_timings_ref_.begin(), 
                                                total_timings_ref_.end()),});
    }
    std::vector<double> ticks(timings.size()/3);
    for (int i = 0; i < ticks.size(); i++){
        ticks[i] = i+1;
    }

    for (int plot = 0; plot < timings.size()/3; plot++){
        std::ofstream file("plotting_data/timings_" + std::to_string(plot) + 
                           ".csv");
        if (file.is_open()){
            file<<"t_solve, t_update, t_total\n";
            for (int i = 0; i < timings[3*plot].size(); i++){
                file<<timings[3*plot][i]<<","<<timings[3*plot+1][i]<<","<<
                    timings[3*plot+2][i]<<"\n";
            }
            file<<getMedian(timings[3*plot])<<","<<
                getMedian(timings[3*plot+1])<<","<<
                getMedian(timings[3*plot+2])<<"\n";

            file.close();
        } 

    }
};

// void Plotter::plotIterCounts(){
//     if (!adaptive_iter_counts_available_ && !ref_iter_counts_available_ &&
//             !naive_iter_counts_available_){
//         return;
//     }

//     plt::figure();
//     if (adaptive_iter_counts_available_){
//         plt::plot(iter_counts_adaptive_, {{"label", "adaptiveNLP"}});
//     }
//     if (ref_iter_counts_available_){
//         plt::plot(iter_counts_ref_, {{"label", "casadi::Opti 1"}});
//     }
//     if (naive_iter_counts_available_){
//         plt::plot(iter_counts_naive_, {{"label", "casadi::Opti 2"}});
//     }
//     plt::legend();

//     plt::save("figures/iter_counts" + file_format_, dpi_);
//     plt::close();
// }

// void Plotter::plotScaling(){
//     plt::figure();

//     if (adaptive_solution_available_){
//         plt::plot({-1000}, {-1000}, {{"color", "royalblue"}, 
//             {"label", "adaptiveNLP"}, {"marker", "o"}});
//         plt::plot(scaling_Ns_adaptive_, scaling_median_times_adaptive_, 
//                 {{"color", "royalblue"}, {"marker", "o"}, 
//                 {"linewidth", "0.4"}});
//         plt::plot(scaling_Ns_adaptive_, scaling_median_total_times_adaptive_, 
//                 {{"color", "royalblue"}, {"marker", "o"}});
//     }
//     if (ref_solution_available_){
//         plt::plot({-1000}, {-1000}, {{"color", "gray"}, 
//             {"label", "casadi::Opti 1"}, {"marker", "o"}});
//         plt::plot(scaling_Ns_ref_, scaling_median_times_ref_, 
//                 {{"color", "gray"}, {"marker", "o"}, {"linewidth", "0.4"}});
//         plt::plot(scaling_Ns_ref_, scaling_median_total_times_ref_, 
//                 {{"color", "gray"}, {"marker", "o"}});
//     }
//     if (naive_solution_available_){
//         plt::plot({-1000}, {-1000}, {{"color", "orange"}, 
//             {"label", "casadi::Opti 2"}, {"marker", "o"}});
//         plt::plot(scaling_Ns_naive_, scaling_median_times_naive_, 
//                 {{"color", "orange"}, {"marker", "o"}, {"linewidth", "0.4"}});
//         plt::plot(scaling_Ns_naive_, scaling_median_total_times_naive_, 
//                 {{"color", "orange"}, {"marker", "o"}});
//     }

//     plt::legend();

//     plt::xlim(scaling_Ns_adaptive_[0]-5, 
//               scaling_Ns_adaptive_[scaling_Ns_adaptive_.size()-1]+5);

//     double min_element = 10000.0;
//     double max_element = -10000.0;
//     std::vector<std::vector<double>> vv = {scaling_median_times_adaptive_, 
//                                            scaling_median_times_naive_, 
//                                            scaling_median_times_ref_, 
//                                            scaling_median_total_times_adaptive_, 
//                                            scaling_median_total_times_naive_, 
//                                            scaling_median_total_times_ref_};
//     for (auto e : vv){
//         for (auto ee : e){
//             min_element = std::min(ee, min_element);
//             max_element = std::max(ee, max_element);
//         }
//     }

//     plt::ylim(min_element-10, max_element+10);

//     plt::xlabel("N");
//     plt::ylabel("time [ms]");

//     plt::save("figures/scaling"+file_format_, 400);
//     plt::close();
// }

void Plotter::writeScaling(){
    assert (adaptive_solution_available_ && naive_solution_available_ && 
        ref_solution_available_);
    
    std::ofstream file("plotting_data/scaling_data.csv");
    if (file.is_open()){
        file<<"N, casadi 1 t_solve, casadi 1 t_total, casadi 2 t_solve, "
            "casadi 2 t_total, adaptiveNLP t_solve, adaptiveNLP t_total\n";
        for (int i = 0; i < scaling_median_times_ref_.size(); i++){
            file<<scaling_Ns_adaptive_[i]<<","<<
                  scaling_median_times_ref_[i]     <<","<<
                    scaling_median_total_times_ref_[i]     <<","<<
                  scaling_median_times_naive_[i]   <<","<<
                    scaling_median_total_times_naive_[i]   <<","<<
                  scaling_median_times_adaptive_[i]<<","<<
                    scaling_median_total_times_adaptive_[i]<<"\n";
        }
        file.close();
    } 
}

void Plotter::addCircle(double x, double y, double r){
    circle_x_.push_back(x);
    circle_y_.push_back(y);
    circle_r_.push_back(r);
    visible_circles_.push_back(true);
    circle_colors_.push_back("r");
}

void Plotter::addCorridor(std::vector<double>& corridor){
    corridors_.push_back(corridor);
}

void Plotter::setCircleVisibility(int ind, bool visibility){
    visible_circles_[ind] = visibility;
}

void Plotter::setCircleColor(int ind, std::string color){
    circle_colors_[ind] = color;
}

void Plotter::formatSolution(std::vector<double>& sol,
                             std::vector<std::vector<double>>& formatted_xx,
                             std::vector<std::vector<double>>& formatted_uu){
    formatted_xx = std::vector<std::vector<double>>(nx_, 
                        std::vector<double>(N_+1, 0.0));
    formatted_uu = std::vector<std::vector<double>>(nu_, 
                        std::vector<double>(N_, 0.0));
    int k = 0;
    for (int i = 0; i < N_+1; i++){
        get_x(&sol[0], k);
        for (int j = 0; j < nx_; j++){
            formatted_xx[j][i] = xk_[j];
        }
        if (k != final_ind_){
            get_u(&sol[0], k);
            for (int j = 0; j < nu_; j++){
                formatted_uu[j][i] = uk_[j];
            }
            k = next_ind_[k].value();
        }
    }
};

void Plotter::printX0s(){
    if (adaptive_solution_available_){
        std::cout<<"Adaptive x0s: "<<std::endl;
        for (auto e : x0s_adaptive_){
            for (auto ee : e){
                std::cout<<ee<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }

    if (ref_solution_available_){
        std::cout<<"Reference x0s: "<<std::endl;
        for (auto e : x0s_ref_){
            for (auto ee : e){
                std::cout<<ee<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }
};

void Plotter::get_x(const double* sol, int k){
	int offset = (k > final_ind_)*nu_;
	for (int i = 0; i < nx_; i++){
		xk_[i] = sol[(nx_ + nu_)*k - offset + free_time_ + i];
	}
}

void Plotter::get_u(const double* sol, int k){
	int offset = (k > final_ind_)*nu_;
	for (int i = 0; i < nu_; i++){
		uk_[i] = sol[(nx_ + nu_)*k + nx_ - offset + free_time_ + i];
	}
}

// void Plotter::plotCircle(int ind){
//     if (visible_circles_[ind]){
//         double x = circle_x_[ind];
//         double y = circle_y_[ind];
//         double r = circle_r_[ind];
//         int N = 100;
//         std::vector<double> xx(N);
//         std::vector<double> yy(N);

//         for (int i = 0; i < N; i++){
//             xx[i] = x + r*std::cos((2*3.141592*i)/(N-1));
//             yy[i] = y + r*std::sin((2*3.141592*i)/(N-1));
//         }
//         plt::plot(xx, yy, {{"color", circle_colors_[ind]}, 
//                            {"label", "obstacle"}});
//     }
// }

// void Plotter::plotCircle(double x, double y, double r, std::string c){
//     int N = 100;
//     std::vector<double> xx(N);
//     std::vector<double> yy(N);

//     for (int i = 0; i < N; i++){
//         xx[i] = x + r*std::cos((2*3.141592*i)/(N-1));
//         yy[i] = y + r*std::sin((2*3.141592*i)/(N-1));
//     }
//     plt::plot(xx, yy, {{"color", c}, {"label", "obstacle"}});
// }

// void Plotter::plotCorridor(int ind){
//     std::vector<double> xx = {corridors_[ind][0], corridors_[ind][1], 
//                               corridors_[ind][1], corridors_[ind][0], 
//                               corridors_[ind][0]};
//     std::vector<double> yy = {corridors_[ind][2], corridors_[ind][2], 
//                               corridors_[ind][3], corridors_[ind][3], 
//                               corridors_[ind][2]};
//     plt::plot(xx, yy, {{"color", "g"}, {"label", "corridor"}});
// }

// void Plotter::plotVehicle(std::vector<double>& corner_points_x, 
//                           std::vector<double>& corner_points_y){
//     std::vector<double> xx = {corner_points_x[0], corner_points_x[1], 
//                               corner_points_x[2], corner_points_x[3], 
//                               corner_points_x[0]};
//     std::vector<double> yy = {corner_points_y[0], corner_points_y[1], 
//                               corner_points_y[2], corner_points_y[3], 
//                               corner_points_y[0]};
//     plt::plot(xx, yy, {{"color", "k"}, {"label", "vehicle"}});
// };

double Plotter::getMedian(std::vector<double>& v){
    std::vector<double> v_sorted = v;
    std::sort(v_sorted.begin(), v_sorted.end());
    
    int size = v_sorted.size();
    if (size & 2 == 1){
        return v_sorted[size/2];
    } else {
        int mid_ind = size/2;
        return (v_sorted[mid_ind-1] + v_sorted[mid_ind])/2.0;
    }
};

std::string Plotter::getMedianText(std::vector<double>& v){
    double median = getMedian(v);
    std::ostringstream oss;
    oss<<std::fixed<<std::setprecision(2)<<std::setw(5)<<
        std::round(median*100.0)/100.0;
    return oss.str();
}; 