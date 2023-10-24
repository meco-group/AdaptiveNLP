#include "core.hpp"
#include "plotter.hpp"
#include "moonlanderExampleHelper.hpp"
#include "errorEstimator.hpp"
#include <chrono>
#include <casadi/casadi.hpp>
#include <fstream>

using namespace std::chrono;
using namespace casadi;

MoonlanderHelper::MoonlanderHelper(int nx_, int nu_){
    nx = nx_; nu = nu_;
    px = std::vector<Polynomial>(nx);
    pu = std::vector<Polynomial>(nu);
};

void MoonlanderHelper::performAdaptiveLoop(BuildingBlocks& blocks, double T, 
                                    int Nmax, double tolerance,
                                    std::vector<double>& solving_times,
                                    std::vector<double>& err_comp_times,
                                    std::vector<double>& updating_times,
                                    std::vector<std::vector<int>>& jac_row,
                                    std::vector<std::vector<int>>& jac_col,
                                    std::vector<std::vector<int>>& hess_row,
                                    std::vector<std::vector<int>>& hess_col,
                                    bool store_sparsities,
                                    int max_iterations,
                                    int print_level){
    
    int nx = blocks.get_nx();
    int nu = blocks.get_nu();

    // create plotter object
    Plotter plotter = Plotter(blocks.getFreeTime(), blocks.get_nx(), 
                              blocks.get_nu());

    // create adaptive NLP object
    AdaptiveNLP adaptiveNLP = AdaptiveNLP(blocks, T, Nmax, 0);
    adaptiveNLP.setIpoptPrintLevel(print_level);
    std::vector<int> nks = {4, 4, 4};

    // get the time-stamps for the collocation grid
    std::vector<double> tau;
    int tt_length = 1;
    for (int i = 0; i < nks.size(); i++){
        tt_length += nks[i]-1;
    }
    std::vector<double> tt(tt_length, 0);
    int tt_ptr = 0;
    for (int k = 0; k < nks.size(); k++){
        tau = casadi::collocation_points(nks[k]-1);
        for (int i = 0; i < tau.size(); i++){
            tt[tt_ptr+i+1] = tt[tt_ptr] + T/nks.size()*tau[i];
        }
        tt_ptr += tau.size();
    }
    adaptiveNLP.initTimeSteps(nks, tt);

    // prepare error estimator
    std::vector<double> error_estimates;
    SX xk = SX::sym("xk", nx, 1);
    SX uk = SX::sym("uk", nu, 1);
    SX p = SX::sym("p", 0, 0);
    Function fd = Function("fd", {xk, uk}, {vertcat(xk(1), -3.8+uk)});
    errorEstimator estimator = errorEstimator(nx, nu, fd);
    std::vector<std::vector<double>> formatted_xx;
    std::vector<std::vector<double>> formatted_uu;
    double max_err = 0.0;
    bool converged = false;
    int iteration_counter = 0;
    int k_offset = 0;
    int curr_interval;
    int new_nb_collocation_pts = 0;
    int new_nb_intervals = 0;
    std::vector<double> tt_local;

    // init guess information
    std::vector<double> init_guess;
    std::vector<std::optional<int>> next_ind_old = adaptiveNLP.getNextInd();
    std::vector<std::optional<int>> nb_g_d_old = adaptiveNLP.getNbSteps();
    std::vector<std::optional<double>> time_from_ind_old = 
        adaptiveNLP.getTimeFromIndex();
    

    // timing containers
    std::vector<double> sol;
    double time;
    if (store_sparsities){
        jac_row = {};
        jac_col = {};
        hess_row = {};
        hess_col = {};
    }

    // adaptive approach
    while (!converged){
        if (store_sparsities){
            jac_row.push_back(adaptiveNLP.getJacRows());
            jac_col.push_back(adaptiveNLP.getJacCols());
            hess_row.push_back(adaptiveNLP.getHessRows());
            hess_col.push_back(adaptiveNLP.getHessCols());
        }

        auto start = high_resolution_clock::now();
        // set initial guess
        if (iteration_counter > 0){
            init_guess = evaluateSolOnNewGrid(sol[0], 
                                              adaptiveNLP.getFinalInd(),
                                              next_ind_old,
                                              nb_g_d_old,
                                              time_from_ind_old, 
                                              adaptiveNLP.getNextInd(),
                                              adaptiveNLP.getNbSteps(),
                                              adaptiveNLP.getTimeFromIndex(),
                                              formatted_xx, formatted_uu);

        } else {
            init_guess = std::vector<double>(1 + nx*(adaptiveNLP.getN()+1) + 
                                             nu*adaptiveNLP.getN(), 0.0);
            init_guess[0] = blocks.get_t_init();
        }
        auto stop = high_resolution_clock::now();
        updating_times.push_back(
            double(duration_cast<microseconds>(stop-start).count())/(1.0e3));

        adaptiveNLP.solveNlp({}, time, init_guess);
        solving_times.push_back(time);
        sol = adaptiveNLP.getSolution();

        // do some post processing
        std::vector<std::optional<int>> next_ind = adaptiveNLP.getNextInd();
        std::vector<std::optional<double>> time_from_ind = 
            adaptiveNLP.getTimeFromIndex();
        plotter.update(adaptiveNLP.getN(), adaptiveNLP.getFinalInd(), 
                    next_ind, time_from_ind);
        plotter.formatSolution(sol, formatted_xx, formatted_uu);

        // compute error estimates
        start = high_resolution_clock::now();
        error_estimates = {};
        int offset = 0;
        k_offset = 0;
        int curr_interval = 0;
        max_err = 0.0;
        while (k_offset != adaptiveNLP.getFinalInd()){

            std::vector<std::vector<double>> xx1(
                adaptiveNLP.getNbSteps(k_offset), std::vector<double>(nx));
            std::vector<std::vector<double>> uu1(
                adaptiveNLP.getNbSteps(k_offset)-1, std::vector<double>(nu));
            for (int i = 0; i < adaptiveNLP.getNbSteps(k_offset); i++){
                for (int j = 0; j < nx; j++){
                    xx1[i][j] = formatted_xx[j][offset + i];
                }
                if (i < adaptiveNLP.getNbSteps(k_offset)-1){
                    for (int j = 0; j < nu; j++){
                        uu1[i][j] = formatted_uu[j][offset + i];
                    }
                }
            }

            offset += adaptiveNLP.getNbSteps(k_offset)-1;
            tau = std::vector<double>(adaptiveNLP.getNbSteps(k_offset));
            tau[0] = 0;
            for (int i = 0; i < adaptiveNLP.getNbSteps(k_offset)-1; i++){
                tau[i+1] = casadi::collocation_points(
                    adaptiveNLP.getNbSteps(k_offset)-1)[i];
            }
            double temp = sol[0]*adaptiveNLP.getIntervalLength(k_offset);

            double err = estimator.getEstimate(
                temp, tau, xx1, uu1);
            max_err = std::max(max_err, err);
            error_estimates.push_back(err);

            int nb_to_traverse = adaptiveNLP.getNbSteps(k_offset)-1; 
            for (int i = 0; i < nb_to_traverse; i++){
                k_offset = adaptiveNLP.getNext(k_offset);
            }
            curr_interval++;

        }
        stop = high_resolution_clock::now();
        time = double(duration_cast<microseconds>(stop-start).count())/(1.0e3);
        err_comp_times.push_back(time);

        // make changes to NLP
        start = high_resolution_clock::now();
        if (max_err <= tolerance){
            stop = high_resolution_clock::now();
            time = double(duration_cast<microseconds>(stop-start).count())/
                    (1.0e3);
            updating_times[updating_times.size()-1] += time;
            converged = true;
        } else {
            k_offset = 0;
            curr_interval = 0;
            while (k_offset != adaptiveNLP.getFinalInd()){
                new_nb_intervals = 1;
                if (error_estimates[curr_interval] > tolerance){
                    new_nb_collocation_pts = 
                        adaptiveNLP.getNbSteps(k_offset) +
                        std::ceil(
                            std::log(error_estimates[curr_interval]/tolerance)/
                            std::log(adaptiveNLP.getNbSteps(k_offset)));
                    
                    // increase number of collocation points if possible
                    if (new_nb_collocation_pts <= 8){
                        tt_local = std::vector<double>(
                                    new_nb_collocation_pts-2);
                        for (int i = 0; i < new_nb_collocation_pts-2; i++){
                            tt_local[i] = 
                                adaptiveNLP.getTimeFromIndex(k_offset) +  
                                adaptiveNLP.getIntervalLength(k_offset)*
                                collocation_points(
                                    new_nb_collocation_pts-1)[i];
                        }
                        adaptiveNLP.changeIntervalDiscretization(
                            k_offset, {new_nb_collocation_pts}, tt_local);
                    
                    // else, split interval
                    } else {
                        new_nb_intervals = std::max(
                            int(std::ceil(new_nb_collocation_pts/4)), 2);
                        double interval_width = 
                            adaptiveNLP.getIntervalLength(k_offset)/
                            new_nb_intervals;
                        tt_local = std::vector<double>(3*new_nb_intervals-1);
                        for (int j = 0; j < new_nb_intervals; j++){
                            for (int i = 0; i < 2 + (j < new_nb_intervals-1); 
                                    i++){
                                tt_local[3*j+i] = 
                                    adaptiveNLP.getTimeFromIndex(k_offset) + 
                                    j*interval_width + 
                                    interval_width*collocation_points(3)[i];
                            }
                        }
                        adaptiveNLP.changeIntervalDiscretization(k_offset, 
                            std::vector<int>(new_nb_intervals, 4), tt_local);
                    }
                }
                curr_interval++;
                for (int j = 0; j < new_nb_intervals; j++){
                    int nb_steps_local = adaptiveNLP.getNbSteps(k_offset);
                    for (int i = 0; i < nb_steps_local-1; i++){
                        k_offset = adaptiveNLP.getNext(k_offset);
                    }
                }
            }
            adaptiveNLP.clearStructuralZeros();

            stop = high_resolution_clock::now();
            time = double(duration_cast<microseconds>(stop-start).count())/
                        (1.0e3);
            updating_times[updating_times.size()-1] += time;

            iteration_counter++;
            if (iteration_counter > max_iterations){
                break;
            }
        }
    }

    std::cout<<"Convergence status:     "<<converged<<endl;
    std::cout<<"number of iterations:   "<<iteration_counter<<endl;

    std::cout<<"solving time:           "<<
        std::accumulate(solving_times.begin(), solving_times.end(), 0.0);
    std::cout<<"   ("<<solving_times<<")"<<endl;
    std::cout<<"error computation time: "<<
        std::accumulate(err_comp_times.begin(), err_comp_times.end(), 0.0);
    std::cout<<"   ("<<err_comp_times<<")"<<endl;
    std::cout<<"updating time:          "<<
        std::accumulate(updating_times.begin(), updating_times.end(), 0.0);
    std::cout<<"   ("<<updating_times<<")"<<endl<<endl<<endl;

    plotter.writeControlsToFile(&sol[0]);
}


void MoonlanderHelper::performCasadiLoop(BuildingBlocks& blocks, double T, 
                                    double tolerance,
                                    std::vector<double>& solving_times,
                                    std::vector<double>& err_comp_times,
                                    std::vector<double>& updating_times,
                                    std::vector<std::vector<int>>& jac_row,
                                    std::vector<std::vector<int>>& jac_col,
                                    std::vector<std::vector<int>>& hess_row,
                                    std::vector<std::vector<int>>& hess_col,
                                    bool store_sparsities,
                                    int max_iterations, 
                                    int print_level){
    int nx = blocks.get_nx();
    int nu = blocks.get_nu();

    std::vector<double> interval_lengths = {1.0/3.0, 1.0/3.0, 1.0/3.0};
    std::vector<int> nb_steps = {4, 4, 4};
    std::vector<double> grid = makeGrid(interval_lengths, nb_steps);
    int N = grid.size()-1;
    std::vector<double> old_grid = grid;
    std::vector<int> nb_steps_old = nb_steps;

    DM xx_DM;
    DM uu_DM;
    double t_double;

    // prepare error estimator
    std::vector<double> error_estimates;
    SX xk = SX::sym("xk", nx, 1);
    SX uk = SX::sym("uk", nu, 1);
    SX p = SX::sym("p", 0, 0);
    Function fd = Function("fd", {xk, uk}, {vertcat(xk(1), -3.8+uk)});
    errorEstimator estimator = errorEstimator(nx, nu, fd);
    std::vector<std::vector<double>> formatted_xx;
    std::vector<std::vector<double>> formatted_uu;
    double max_err = 0.0;
    bool converged = false;
    int iteration_counter = 0;
    int k_offset = 0;
    int curr_interval;
    int new_nb_collocation_pts = 0;
    int new_nb_intervals = 0;
    std::vector<double> tau;

    double time;
    
    if (store_sparsities){
        jac_row = {};
        jac_col = {};
        hess_row = {};
        hess_col = {};
    }


    while (!converged){
        if (store_sparsities){
            jac_row.push_back({});
            jac_col.push_back({});
            hess_row.push_back({});
            hess_col.push_back({});
        }
        auto start = high_resolution_clock::now();
        Opti opti = Opti();
        MX xx = opti.variable(nx, N+1);
        MX uu = opti.variable(nu, N);
        MX t = opti.variable(1, 1);

        // objective
        MX obj = blocks.eval_Phi_0({xx(Slice(),0), t, 0})[0] + 
                blocks.eval_Phi_f({xx(Slice(),N), t, 0})[0];
        for (int k = 0; k < N; k++){
            obj += blocks.eval_phi({xx(Slice(), k), uu(Slice(), k), 
                            (grid[k+1]-grid[k]), t, 0})[0];
        }
        opti.minimize(obj);

        // constraints
        opti.subject_to((blocks.get_g0_lb() <= 
                        blocks.eval_g0({xx(Slice(),0), t, 0})[0]) <= 
                        blocks.get_g0_ub());
        for (int k = 0; k < N; k++){
            opti.subject_to((blocks.get_g_fixed_lb() <= 
                blocks.eval_g_fixed({xx(Slice(), k), uu(Slice(), k), 0})[0]) <= 
                blocks.get_g_fixed_ub());
        }

        opti.subject_to((blocks.get_gT_lb() <= 
                        blocks.eval_gT({xx(Slice(),N), 0})[0]) 
                        <= blocks.get_gT_ub());

        int k_offset = 0;
        std::vector<double> dts;
        for (int k = 0; k < nb_steps.size(); k++){
            MX local_xx = MX::zeros(nx, nb_steps[k]);
            MX local_uu = MX::zeros(nu, nb_steps[k]-1);
            for (int j = 0; j < nb_steps[k]-1; j++){
                local_xx(Slice(), j) = xx(Slice(), k_offset+j);
                local_uu(Slice(), j) = uu(Slice(), k_offset+j);
            }
            local_xx(Slice(), nb_steps[k]-1) = xx(Slice(), 
                     k_offset+nb_steps[k]-1);
            
            dts = std::vector<double>(nb_steps[k]-1);
            for (int i = 0; i < dts.size(); i++){
                dts[i] = grid[k_offset+i+1] - grid[k_offset+i];
            }

            opti.subject_to(
                blocks.eval_g_disc
                    ({local_xx, local_uu, dts, t, 0}, nb_steps[k])[0] == 0);

            k_offset += nb_steps[k]-1;
        }

        // initial guess
        if (iteration_counter > 0){
            std::vector<std::vector<double>> init_guess = 
                evaluateSolOnNewGrid(old_grid, nb_steps_old, formatted_xx, 
                                     formatted_uu, grid);
            for (int k = 0; k < N; k++){
                opti.set_initial(xx(0, k), init_guess[0][k]);
                opti.set_initial(xx(1, k), init_guess[1][k]);
                if (k < N-1){
                    opti.set_initial(uu(0, k), init_guess[2][k+1]);
                }
            }
            opti.set_initial(t, t_double);
        } else {
            for (int i = 0; i < N+1; i++){ 
                opti.set_initial(xx(Slice(), i), 
                                 blocks.eval_x_init({DM(grid[i])})[0]);
            }
            for (int i = 0; i < N; i++){ 
                opti.set_initial(uu(Slice(), i), 
                                 blocks.eval_u_init({DM(grid[i])})[0]);
            }
            opti.set_initial(t, 5);
        }

        // set parameters
        opti.solver("ipopt", {{"print_time", false}, {"record_time", true}}, {{"print_level", print_level}});
        auto stop = high_resolution_clock::now();
        updating_times.push_back(double(
            duration_cast<microseconds>(stop-start).count())/(1.0e3)); 

        // solve problem
        OptiSol sol_casadi = opti.solve();
        solving_times.push_back(1000*double(sol_casadi.stats()["t_wall_total"]));

        if (store_sparsities){
            getSparsities(opti, jac_row[iteration_counter], 
                        jac_col[iteration_counter], hess_row[iteration_counter],
                        hess_col[iteration_counter]);
        }

        xx_DM = sol_casadi.value(xx);
        uu_DM = sol_casadi.value(uu);
        t_double = double(sol_casadi.value(t));
        formatted_xx = std::vector<std::vector<double>>(nx, 
                            std::vector<double>(N+1));
        formatted_uu = std::vector<std::vector<double>>(nu, 
                            std::vector<double>(N));
        for (int i = 0; i < nx; i++){
            for (int j = 0; j < N+1; j++){
                formatted_xx[i][j] = double(xx_DM(i,j));
            }
        }
        for (int i = 0; i < nu; i++){
            for (int j = 0; j < N; j++){
                formatted_uu[i][j] = double(uu_DM(i,j));
            }
        }

        // compute error estimates
        start = high_resolution_clock::now();
        error_estimates = {};
        int offset = 0;
        max_err = 0.0;
        for (int curr_interval = 0; curr_interval < nb_steps.size(); 
                curr_interval++){
            int nb_steps_local = nb_steps[curr_interval];
            
            std::vector<std::vector<double>> xx1(
                nb_steps_local, std::vector<double>(nx));
            std::vector<std::vector<double>> uu1(
                nb_steps_local-1, std::vector<double>(nu));
            for (int i = 0; i < nb_steps_local; i++){
                for (int j = 0; j < nx; j++){
                    xx1[i][j] = formatted_xx[j][offset + i];
                }
                if (i < nb_steps_local-1){
                    for (int j = 0; j < nu; j++){
                        uu1[i][j] = formatted_uu[j][offset + i];
                    }
                }
            }

            offset += nb_steps_local-1;
            tau = std::vector<double>(nb_steps_local);
            tau[0] = 0;
            for (int i = 0; i < nb_steps_local-1; i++){
                tau[i+1] = casadi::collocation_points(
                    nb_steps_local-1)[i];
            }
            double temp = t_double*interval_lengths[curr_interval];

            double err = estimator.getEstimate(
                temp, tau, xx1, uu1);
            max_err = std::max(max_err, err);
            error_estimates.push_back(err);
        }
        stop = high_resolution_clock::now();
        time = double(duration_cast<microseconds>(stop-start).count())/(1.0e3);
        err_comp_times.push_back(time);

        // make changes to NLP
        old_grid = grid;
        nb_steps_old = nb_steps;
        start = high_resolution_clock::now();
        if (max_err <= tolerance){
            stop = high_resolution_clock::now();
            time = double(duration_cast<microseconds>(stop-start).count())/
                    (1.0e3);
            updating_times[updating_times.size()-1] += time;
            converged = true;
        } else {
            std::vector<int> nb_steps_new = {};
            std::vector<double> interval_lengths_new = {};

            for (int curr_interval = 0; curr_interval < nb_steps.size();
                    curr_interval++){
                new_nb_intervals = 1;
                if (error_estimates[curr_interval] > tolerance){
                    new_nb_collocation_pts = 
                        nb_steps[curr_interval] +
                        std::ceil(
                            std::log(error_estimates[curr_interval]/tolerance)/
                            std::log(nb_steps[curr_interval]));
                   
                    // increase number of collocation points if possible
                    if (new_nb_collocation_pts <= 8){                        
                        nb_steps_new.push_back(new_nb_collocation_pts);
                        interval_lengths_new.push_back(
                            interval_lengths[curr_interval]);
                    
                    // else, split interval
                    } else {
                        new_nb_intervals = std::max(
                            int(std::ceil(new_nb_collocation_pts/4)), 2);
                        double interval_width = 
                            interval_lengths[curr_interval]/new_nb_intervals;

                        for (int j = 0; j < new_nb_intervals; j++){
                            nb_steps_new.push_back(4);
                            interval_lengths_new.push_back(interval_width);
                        }
                    }
                } else {
                    nb_steps_new.push_back(nb_steps[curr_interval]);
                    interval_lengths_new.push_back(
                        interval_lengths[curr_interval]);
                }
            }

            nb_steps = nb_steps_new;
            interval_lengths = interval_lengths_new;
            grid = makeGrid(interval_lengths, nb_steps);
            N = grid.size()-1;

            stop = high_resolution_clock::now();
            time = double(duration_cast<microseconds>(stop-start).count())/
                    (1.0e3);
            updating_times[updating_times.size()-1] += time;

            iteration_counter++;
            if (iteration_counter > max_iterations){
                break;
            }
        }

    }

    std::cout<<"Convergence status:     "<<converged<<endl;
    std::cout<<"number of iterations:   "<<iteration_counter<<endl;

    std::cout<<"solving time:           "<<
        std::accumulate(solving_times.begin(), solving_times.end(), 0.0);
    std::cout<<"   ("<<solving_times<<")"<<endl;
    std::cout<<"error computation time: "<<
        std::accumulate(err_comp_times.begin(), err_comp_times.end(), 0.0);
    std::cout<<"   ("<<err_comp_times<<")"<<endl;
    std::cout<<"updating time:          "<<
        std::accumulate(updating_times.begin(), updating_times.end(), 0.0);
    std::cout<<"   ("<<updating_times<<")"<<endl;
}

void MoonlanderHelper::writeToFile(std::vector<double>& solving_times_adaptive,
                                std::vector<double>& err_comp_times_adaptive,
                                std::vector<double>& updating_times_adaptive,
                                std::vector<double>& solving_times_casadi,
                                std::vector<double>& err_comp_times_casadi,
                                std::vector<double>& updating_times_casadi){
    std::ofstream file("../../examples/plotting_data/adaptive_gridding_timings.csv");
    if (file.is_open()){
        file << "t_solve_adaptive, t_err_adative, t_update_adaptive, "
                "t_solve_casadi, t_err_casadi, t_update_casadi\n";
        for (int i = 0; i < solving_times_adaptive.size(); i++){
            file<<solving_times_adaptive[i]<<", ";
            file<<err_comp_times_adaptive[i]<<", ";
            file<<updating_times_adaptive[i]<<", ";
            file<<solving_times_casadi[i]<<", ";
            file<<err_comp_times_casadi[i]<<", ";
            file<<updating_times_casadi[i]<<"\n";
        }
        file.close();
    }
}

void MoonlanderHelper::writeToFile(
                         std::vector<std::vector<int>>& jac_rows_adaptive,
                         std::vector<std::vector<int>>& jac_cols_adaptive,
                         std::vector<std::vector<int>>& hess_rows_adaptive,
                         std::vector<std::vector<int>>& hess_cols_adaptive,
                         std::vector<std::vector<int>>& jac_rows_casadi,
                         std::vector<std::vector<int>>& jac_cols_casadi,
                         std::vector<std::vector<int>>& hess_rows_casadi,
                         std::vector<std::vector<int>>& hess_cols_casadi){
    std::ofstream file1("../../examples/plotting_data/adaptive_gridding_sparsities_jac.csv");
    if (file1.is_open()){
        file1 << "jac_rows_adaptive, jac_cols_adative, jac_rows_casadi, "
                 "jac_cols_casadi\n";
        for (int i = 0; i < jac_rows_adaptive.size(); i++){
            for (int j = 0; j < jac_rows_adaptive[i].size(); j++){
                file1<<jac_rows_adaptive[i][j]<<", ";
                file1<<jac_cols_adaptive[i][j]<<", ";
                file1<<jac_rows_casadi[i][j]<<", ";
                file1<<jac_cols_casadi[i][j]<<"\n";
            }
            file1<<"-1,-1,-1,-1\n";
        }
        file1.close();
    }

    std::ofstream file2("../../examples/plotting_data/adaptive_gridding_sparsities_hess.csv");
    if (file2.is_open()){
        file2 << "hess_rows_adaptive, hess_cols_adaptive, "
                "hess_rows_casadi, hess_cols_casadi\n";
        for (int i = 0; i < hess_rows_adaptive.size(); i++){
            for (int j = 0; j < hess_rows_adaptive[i].size(); j++){
                file2<<hess_rows_adaptive[i][j]<<", ";
                file2<<hess_cols_adaptive[i][j]<<", ";
                file2<<hess_rows_casadi[i][j]<<", ";
                file2<<hess_cols_casadi[i][j]<<"\n";
            }
            file2<<"-1,-1,-1,-1\n";
        }
        file2.close();
    }
}

std::vector<double> MoonlanderHelper::makeGrid(
                                        std::vector<double>& interval_lengths, 
                                        std::vector<int>& nb_steps){
    assert (interval_lengths.size() == nb_steps.size());
    int total_nb_steps = 1;
    for (int k = 0; k < interval_lengths.size(); k++){
        total_nb_steps += nb_steps[k]-1;
    }
    
    std::vector<double> grid(total_nb_steps, 0.0);

    int grid_ptr = 0;
    std::vector<double> tau;
    for (int k = 0; k < interval_lengths.size(); k++){
        tau = collocation_points(nb_steps[k]-1);
        for (int i = 0; i < nb_steps[k]-1; i++){
            grid[grid_ptr+i+1] = grid[grid_ptr] + interval_lengths[k]*tau[i];
        }
        grid_ptr += nb_steps[k]-1;
    }

    return grid;
}

std::vector<std::vector<double>> MoonlanderHelper::evaluateSolOnNewGrid(
                                std::vector<double>& old_grid,
                                std::vector<int>& nb_steps_old,
                                std::vector<std::vector<double>>& formatted_xx,
                                std::vector<std::vector<double>>& formatted_uu,
                                std::vector<double>& new_grid){
    std::vector<std::vector<double>> xx_local;
    std::vector<std::vector<double>> uu_local;
    std::vector<double> tau_local;

    std::vector<std::vector<double>> res(nx+nu, 
        std::vector<double>(new_grid.size()));

    int new_k_ptr = 0;

    int offset = 0;
    double interval_edge = 0.0;
    for (int k = 0; k < nb_steps_old.size(); k++){
        xx_local = std::vector<std::vector<double>>(nx,
                        std::vector<double>(nb_steps_old[k]));
        uu_local = std::vector<std::vector<double>>(nu,
                        std::vector<double>(nb_steps_old[k]-1));
        for (int i = 0; i < nb_steps_old[k]; i++){
            for (int j = 0; j < nx; j++){
                xx_local[j][i] = formatted_xx[j][offset+i];
            }
            if (i < nb_steps_old[k]-1){
                for (int j = 0; j < nu; j++){
                    uu_local[j][i] = formatted_uu[j][offset+i];
                }
            }
        }
        
        tau_local = std::vector<double>(nb_steps_old[k]-1);
        for (int i = 0; i < tau_local.size(); i++){
            tau_local[i] = old_grid[offset + i];
        }
        interval_edge = tau_local[tau_local.size()-1];
        offset += nb_steps_old[k]-1;

        constructInterpolants(tau_local, xx_local, uu_local);

        while (new_grid[new_k_ptr] <= interval_edge){
            for (int i = 0; i < nx; i++){
                res[i][new_k_ptr] = px[i](new_grid[new_k_ptr]);
            }
            for (int i = 0; i < nu; i++){
                res[nx+i][new_k_ptr] = pu[i](new_grid[new_k_ptr]);
            }
            new_k_ptr++;
        }
    }

    return res;
};

std::vector<double> MoonlanderHelper::evaluateSolOnNewGrid(double t_old,
                        int final_ind,
                        std::vector<std::optional<int>>& next_ind_old,
                        std::vector<std::optional<int>>& nb_g_d_old,
                        std::vector<std::optional<double>>& time_from_ind_old,
                        std::vector<std::optional<int>> next_ind_new,
                        std::vector<std::optional<int>> nb_g_d_new,
                        std::vector<std::optional<double>> time_from_ind_new,
                        std::vector<std::vector<double>>& formatted_xx,
                        std::vector<std::vector<double>>& formatted_uu){
   
    std::vector<std::vector<double>> xx_local;
    std::vector<std::vector<double>> uu_local;
    std::vector<double> tau_local;

    int N_old = 0;
    int ptr = 0;
    while(next_ind_old[ptr].has_value()){
        N_old++;
        ptr = next_ind_old[ptr].value();
    }

    int N_new = 0;
    ptr = 0;
    while(next_ind_new[ptr].has_value()){
        N_new++;
        ptr = next_ind_new[ptr].value();
    }

    std::vector<double> res(1 + nx*(N_new+1) + nu*N_new);
    res[0] = t_old;

    int new_k_ptr = 0;

    int offset = 0;
    double interval_edge = 0.0;
    int k = 0;
    while (next_ind_old[k].has_value()){
        xx_local = std::vector<std::vector<double>>(nx,
                        std::vector<double>(nb_g_d_old[k].value()));
        uu_local = std::vector<std::vector<double>>(nu,
                        std::vector<double>(nb_g_d_old[k].value()-1));
        for (int i = 0; i < nb_g_d_old[k].value(); i++){
            for (int j = 0; j < nx; j++){
                xx_local[j][i] = formatted_xx[j][offset+i];
            }
            if (i < nb_g_d_old[k].value()-1){
                for (int j = 0; j < nu; j++){
                    uu_local[j][i] = formatted_uu[j][offset+i];
                }
            }
        }
        
        tau_local = std::vector<double>(nb_g_d_old[k].value()-1);
        int temp = k;
        for (int i = 0; i < nb_g_d_old[k].value()-1; i++){
            tau_local[i] = time_from_ind_old[temp].value();
            temp = next_ind_old[temp].value();
        }
        interval_edge = tau_local[tau_local.size()-1];
        offset += nb_g_d_old[k].value()-1;

        constructInterpolants(tau_local, xx_local, uu_local);

        while (time_from_ind_new[new_k_ptr].value() <= interval_edge){
            for (int i = 0; i < nx; i++){
                res[1+(nx+nu)*new_k_ptr-nu*(new_k_ptr > final_ind)+i] = 
                    px[i](time_from_ind_new[new_k_ptr].value());
            }
            for (int i = 0; i < nu; i++){
                res[1+(nx+nu)*new_k_ptr-nu*(new_k_ptr > final_ind)+nx+i] = 
                    pu[i](time_from_ind_new[new_k_ptr].value());
            }
            new_k_ptr = next_ind_new[new_k_ptr].value();
        }

        int nb_steps = nb_g_d_old[k].value();
        for (int i = 0; i < nb_steps-1; i++){
            k = next_ind_old[k].value();
        }
    }

    return res;
};

void MoonlanderHelper::constructInterpolants(std::vector<double>& old_grid,
                            std::vector<std::vector<double>>& formatted_xx,
                            std::vector<std::vector<double>>& formatted_uu){
    std::vector<double> etau = old_grid; 
    if (etau[0] > 0){etau.insert(etau.begin(), 0.);}
    Polynomial pLagrange;
    // for every state
    for (int k = 0; k < nx; k++){
        // Construct interplating polynomial by first constructing lagrange 
        // polynomials
        px[k] = 0;
        // for every lagrange polynomial
        for (int j = 0; j < etau.size(); j++){
            pLagrange = 1;
            // construct single Lagrange polynomial
            for (int r = 0; r < etau.size(); r++){
                if (r != j){
                    pLagrange *= Polynomial(-etau[r],1)/(etau[j]-etau[r]);
                }
            }
            // Add weighted lagrange polynomial to the interpolant
            pLagrange *= Polynomial(formatted_xx[k][j]);
            px[k] += pLagrange;
        }
        if (px[k].degree() < 0){
            px[k] = 0;
        }
    }

    // for every control
    for (int k = 0; k < nu; k++){
        // Construct interplating polynomial by first constructing lagrange 
        // polynomials
        pu[k] = 0;
        // for every lagrange polynomial
        for (int j = 1; j < etau.size(); j++){
            pLagrange = 1;
            // construct single Lagrange polynomial
            for (int r = 1; r < etau.size(); r++){
                if (r != j){
                    pLagrange *= Polynomial(-etau[r],1)/(etau[j]-etau[r]);
                }
            }
            // Add weighted lagrange polynomial to the interpolant
            pLagrange *= Polynomial(formatted_uu[k][j-1]);
            pu[k] += pLagrange;
        }
        if (pu[k].degree() < 0){
            pu[k] = 0;
        }
    }
};

void MoonlanderHelper::getSparsities(Opti opti, std::vector<int>& jac_row, 
                                     std::vector<int>& jac_col, 
                                     std::vector<int>& hess_row,
                                     std::vector<int>& hess_col){
    MX J = jacobian(opti.g(), opti.x());
    Sparsity sp = J.sparsity();
    // // Function J = opti.advanced().get_solver().get_function("nlp_jac_g");
    // Sparsity sp = J.sparsity_out(1);

    for (int i = 0; i < sp.size1(); i++){
        for (int j = 0; j < sp.size2(); j++){
            if (sp.has_nz(i,j)){
                jac_row.push_back(i);
                jac_col.push_back(j);
            }
        }
    }

    MX H = hessian(opti.f() + mtimes(transpose(opti.lam_g()), opti.g()), opti.x());
    sp = H.sparsity();
    // Function H = opti.advanced().get_solver().get_function("nlp_hess_l");
    // sp = H.sparsity_out(0);

    for (int i = 0; i < sp.size1(); i++){
        for (int j = i; j < sp.size2(); j++){
            if (sp.has_nz(i,j)){
                hess_row.push_back(i);
                hess_col.push_back(j);
            }
        }
    }
}