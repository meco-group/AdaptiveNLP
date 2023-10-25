#include "core.hpp"
#include "makeBuildingBlocks.hpp"
#include "plotter.hpp"
#include "errorEstimator.hpp"
#include "moonlanderExampleHelper.hpp"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main(){
    //////////////////////////////////////
    // general definition of parameters //
    //////////////////////////////////////

    double T = 1.0;
    int Nmax = 100;
    double tolerance = 1.0e-6;

    /////////////////////////////
    // create object instances //
    /////////////////////////////

    BuildingBlocks blocks = makeBlocksMoonlander();
    cout<<"BuildingBlocks created!"<<endl;
    
    MoonlanderHelper helper = MoonlanderHelper(blocks.get_nx(), 
                                               blocks.get_nu());
    
    std::string file_path = __FILE__;
    std::string dir_path = file_path.substr(0, file_path.rfind("/"));
    helper.setOutputFolder(dir_path + "/plotting_data/");
    
    int max_nb_refinements = 10;
    int print_level = 0;

    // containers for timings
    std::vector<double> t_solve_adaptive = {};
    std::vector<double> t_err_adaptive = {};
    std::vector<double> t_update_adaptive = {};

    std::vector<double> t_solve_casadi = {};
    std::vector<double> t_err_casadi = {};
    std::vector<double> t_update_casadi = {};

    // containers for sparsities
    std::vector<std::vector<int>> jac_rows_adaptive = {};
    std::vector<std::vector<int>> jac_cols_adaptive = {};
    std::vector<std::vector<int>> hess_rows_adaptive = {};
    std::vector<std::vector<int>> hess_cols_adaptive = {};

    std::vector<std::vector<int>> jac_rows_casadi = {};
    std::vector<std::vector<int>> jac_cols_casadi = {};
    std::vector<std::vector<int>> hess_rows_casadi = {};
    std::vector<std::vector<int>> hess_cols_casadi = {};

    int nb_runs = 100;

    cout<<endl<<endl<<"ADAPTIVE"<<endl;
    for (int counter = 0; counter < nb_runs; counter++){
        helper.performAdaptiveLoop(blocks, T, Nmax, tolerance, 
                                   t_solve_adaptive, t_err_adaptive, 
                                   t_update_adaptive, jac_rows_adaptive, 
                                   jac_cols_adaptive, hess_rows_adaptive,
                                   hess_cols_adaptive, counter==0, 
                                   max_nb_refinements, print_level);
    }
    cout<<endl<<endl<<"CASADI"<<endl;
    for (int counter = 0; counter < nb_runs; counter++){
        helper.performCasadiLoop(blocks, T, tolerance, t_solve_casadi, 
                                t_err_casadi, t_update_casadi, 
                                jac_rows_casadi, jac_cols_casadi, 
                                hess_rows_casadi, hess_cols_casadi,
                                counter==0, max_nb_refinements, print_level);
    }

    helper.writeToFile(t_solve_adaptive, t_err_adaptive, t_update_adaptive,
                       t_solve_casadi, t_err_casadi, t_update_casadi);

    helper.writeToFile(jac_rows_adaptive, jac_cols_adaptive, 
                       hess_rows_adaptive, hess_cols_adaptive, jac_rows_casadi,
                       jac_cols_casadi, hess_rows_casadi, hess_cols_casadi);

    return 0;
}