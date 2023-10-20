#ifndef __MANY_OBSTACLES_HELPER__
#define __MANY_OBSTACLES_HELPER__
#include "adaptiveNLP.hpp"
#include "plotter.hpp"

class ManyObstaclesHelper{
    public:
        ManyObstaclesHelper(AdaptiveNLP& nlp);

        void updateNLPVars();

        std::vector<int> getCollidingPoints(double pos_x, double pos_y, 
                                            double r,
                                            std::vector<double>& sol, 
                                            double margin);
                                        
        bool checkFeasible(std::vector<double>& sol, std::vector<double>& xx,
                           std::vector<double>& yy, std::vector<double>& rr);
        
        bool processObstacleConstraints(std::vector<double>& sol, 
                                        std::vector<double> xx, 
                                        std::vector<double> yy,
                                        std::vector<double> rr);

        bool processObstacleConstraintsCasadi(casadi::Opti& opti, 
                                        MX& xx, MX& uu, BuildingBlocks& blocks,
                                        std::vector<std::vector<double>>& sol, 
                                        std::vector<double> x, 
                                        std::vector<double> y,
                                        std::vector<double> r);

        double performAdaptiveLoop(Plotter& plotter, 
            std::vector<double>& x0, std::vector<double>& xf, int nx, int nu, 
            std::vector<double>& obstacle_x, std::vector<double>& obstacle_y, 
            std::vector<double>& obstacle_r, int max_nb_iterations, 
            std::vector<double>& timings_adaptive, 
            std::vector<std::vector<std::vector<double>>>& 
                formatted_solutions,
            std::vector<int>& nb_constraints, double& total_time,
            int& iteration_count);

        double solveCompleteNLPCasadi(BuildingBlocks& blocks, Plotter& plotter, 
            std::vector<double> x0, std::vector<double> xf, int nx, int nu, 
            int nb_steps, double dt, int nb_intervals, 
            std::vector<double>& obstacle_x, std::vector<double>& obstacle_y, 
            std::vector<double>& obstacle_r, 
            std::vector<std::vector<std::vector<double>>>& 
                formatted_solutions, double& total_time);

        double performCasadiLoopReformulation(BuildingBlocks& blocks, 
            Plotter& plotter, std::vector<double> x0, std::vector<double> xf, 
            int nx, int nu, int nb_steps, double dt, int nb_intervals, 
            std::vector<double>& obstacle_x, std::vector<double>& obstacle_y, 
            std::vector<double>& obstacle_r, int max_nb_iterations,
            std::vector<double>& timings_casadi_2,
            std::vector<std::vector<std::vector<double>>>& 
                formatted_solutions,
            double& total_time, int& iteration_count);
    
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

        std::vector<double> xk_;
        std::vector<double> uk_;
        std::vector<int> all_;
};

#endif