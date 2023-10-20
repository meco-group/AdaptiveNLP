#ifndef __MOONLANDER_HELPER__
#define __MOONLANDER_HELPER__
#include "buildingBlocks.hpp"

class MoonlanderHelper {
    public:
        MoonlanderHelper(int nx_, int nu_);

        void performAdaptiveLoop(BuildingBlocks& blocks, double T, int Nmax, 
                                 double tolerance, 
                                 std::vector<double>& solving_times,
                                 std::vector<double>& err_comp_times,
                                 std::vector<double>& updating_times,
                                 std::vector<std::vector<int>>& jac_row,
                                 std::vector<std::vector<int>>& jac_col,
                                 std::vector<std::vector<int>>& hess_row,
                                 std::vector<std::vector<int>>& hess_col,
                                 bool store_sparsities,
                                 int max_iterations=10,
                                 int print_level=5);

        void performCasadiLoop(BuildingBlocks& blocks, double T, 
                               double tolerance, 
                               std::vector<double>& solving_times,
                               std::vector<double>& err_comp_times,
                               std::vector<double>& updating_times,
                               std::vector<std::vector<int>>& jac_row,
                               std::vector<std::vector<int>>& jac_col,
                               std::vector<std::vector<int>>& hess_row,
                               std::vector<std::vector<int>>& hess_col,
                               bool store_sparsities,
                               int max_iterations=10,
                               int print_level=5);

        void writeToFile(std::vector<double>& solving_times_adaptive,
                         std::vector<double>& err_comp_times_adaptive,
                         std::vector<double>& updating_times_adaptive,
                         std::vector<double>& solving_times_casadi,
                         std::vector<double>& err_comp_times_casadi,
                         std::vector<double>& updating_times_casadi);

        void writeToFile(std::vector<std::vector<int>>& jac_rows_adaptive,
                         std::vector<std::vector<int>>& jac_cols_adaptive,
                         std::vector<std::vector<int>>& hess_rows_adaptive,
                         std::vector<std::vector<int>>& hess_cols_adaptive,
                         std::vector<std::vector<int>>& jac_rows_casadi,
                         std::vector<std::vector<int>>& jac_cols_casadi,
                         std::vector<std::vector<int>>& hess_rows_casadi,
                         std::vector<std::vector<int>>& hess_cols_casadi);


    private:
        std::vector<double> makeGrid(std::vector<double>& interval_lengths, 
                                     std::vector<int>& nb_steps);

        std::vector<std::vector<double>> evaluateSolOnNewGrid(
                                std::vector<double>& old_grid,
                                std::vector<int>& nb_steps_old,
                                std::vector<std::vector<double>>& formatted_xx,
                                std::vector<std::vector<double>>& formatted_uu,
                                std::vector<double>& new_grid);

        std::vector<double> evaluateSolOnNewGrid(double t_old, int final_ind,
                        std::vector<std::optional<int>>& next_ind_old,
                        std::vector<std::optional<int>>& nb_g_d_old,
                        std::vector<std::optional<double>>& time_from_ind_old,
                        std::vector<std::optional<int>> next_ind_new,
                        std::vector<std::optional<int>> nb_g_d_new,
                        std::vector<std::optional<double>> time_from_ind_new,
                        std::vector<std::vector<double>>& formatted_xx,
                        std::vector<std::vector<double>>& formatted_uu);
        
        void constructInterpolants(std::vector<double>& old_grid,
                            std::vector<std::vector<double>>& formatted_xx,
                            std::vector<std::vector<double>>& formatted_uu);

        void getSparsities(Opti opti, std::vector<int>& jac_row, 
                           std::vector<int>& jac_col, 
                           std::vector<int>& hess_row,
                           std::vector<int>& hess_col);

        int nx;
        int nu;
        std::vector<Polynomial> px;
        std::vector<Polynomial> pu;
};

#endif