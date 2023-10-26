#ifndef __CORRIDOR_EXAMPLE_HELPERS_
#define __CORRIDOR_EXAMPLE_HELPERS_
#include <vector>
#include "core.hpp"
#include "plotter.hpp"
#include <optional>

// Helper object to perform the MPC-loop for the different cases
class AdaptiveCorridorHelper{
    public:
        AdaptiveCorridorHelper(AdaptiveNLP& nlp, double view_radius,
                               double low_speed_min);

        void updateNLPVars();

        void getCollidingPoints(double pos_x, double pos_y, double r, 
                                std::vector<double>& sol,
                                std::vector<bool>& constraint_added,
                                std::vector<int>& to_add,
                                std::vector<int>& to_remove);

        std::vector<int> getLeavingPoints(std::vector<double>& corridor, 
                                            std::vector<double>& sol);

        std::vector<int> getLeavingPoints(std::vector<double>& corridor,
                            std::vector<std::vector<double>>& formatted_xx);

        bool processObstacleConstraints(std::vector<double>& x0,
                                          std::vector<double>& sol, 
                                          double pos_x, double pos_y, double r,
                                          int c_idx, int instance_idx);
        
        bool processObstacleConstraintsCoarse(std::vector<double>& x0, 
                                              double pos_x, double pos_y, 
                                              double r, int c_idx,  
                                              int instance_idx, bool seen);

        bool processLowSpeedZone(std::vector<double>& x0,
                                    std::vector<double>& sol, 
                                    double pos_x, double pos_y, double r,
                                    int c_idx, int instance_idx);

        bool processLowSpeedZoneCoarse(std::vector<double>& x0,
                                       double pos_x, double pos_y, double r,
                                       int c_idx, int instance_idx, bool seen);

        std::vector<double> shiftInitialGuess(std::vector<double>& sol);

        void appendInitGuess(std::vector<double>& sol, 
                             std::vector<double>& init_guess, int nb_steps);

        std::vector<std::vector<double>> performAdaptiveLoop(
            AdaptiveNLP& adaptiveNLP, Plotter& plotter, 
            std::vector<double> x0, std::vector<double>& corridor_1, 
            std::vector<double>& corridor_2, int nx, int nu, int nb_steps, 
            int nb_iterations, double dt, int nb_intervals, 
            std::vector<double>& obstacle_x, std::vector<double>& obstacle_y, 
            std::vector<double>& obstacle_r, 
            std::vector<double>& low_speed_pos_x, 
            std::vector<double>& low_speed_pos_y, 
            std::vector<double>& low_speed_r);

        std::vector<std::vector<double>> performCasadiLoop(
            BuildingBlocks& blocks, Plotter& plotter, 
            std::vector<double> x0, std::vector<double>& corridor_1, 
            std::vector<double>& corridor_2, int nx, int nu, int nb_steps, 
            int nb_iterations, double dt, double T, int nb_intervals, 
            std::vector<double>& obstacle_x, std::vector<double>& obstacle_y, 
            std::vector<double>& obstacle_r, 
            std::vector<double>& low_speed_pos_x, 
            std::vector<double>& low_speed_pos_y, 
            std::vector<double>& low_speed_r);

        std::vector<std::vector<double>> performCasadiLoopReformulation(
            BuildingBlocks& blocks, Plotter& plotter, 
            std::vector<double> x0, std::vector<double>& corridor_1, 
            std::vector<double>& corridor_2, int nx, int nu, int nb_steps, 
            int nb_iterations, double dt, double T, int nb_intervals, 
            std::vector<double>& obstacle_x, std::vector<double>& obstacle_y, 
            std::vector<double>& obstacle_r, 
            std::vector<double>& low_speed_pos_x, 
            std::vector<double>& low_speed_pos_y, 
            std::vector<double>& low_speed_r);

    private:
        void get_x(std::vector<double>& sol, int k);

        void get_u(std::vector<double>& sol, int k);

        AdaptiveNLP* const nlp_;

        int nx_;
        int nu_;
        int N_;
        std::vector<std::optional<int>> next_index_;
        std::vector<std::optional<double>> time_from_index_;
        int final_ind_;
        bool free_time_;

        double view_radius_;
        double low_speed_min_;

        std::vector<double> xk_;
        std::vector<double> uk_;
        std::vector<int> all_;

        std::vector<std::vector<bool>> obs_constraint_added_;
        std::vector<std::vector<bool>> low_constraint_added_;
};

#endif