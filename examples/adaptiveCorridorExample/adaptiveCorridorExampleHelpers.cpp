#include "adaptiveCorridorExampleHelpers.hpp"
#include "core.hpp"
#include <vector>
#include <cmath>
#include <casadi/casadi.hpp>
#include <chrono>

using namespace std::chrono;

AdaptiveCorridorHelper::AdaptiveCorridorHelper(AdaptiveNLP& nlp,
                                               double view_radius,
                                               double low_speed_min)
                                               :nlp_(&nlp){
    view_radius_ = view_radius;
    low_speed_min_ = low_speed_min;
    updateNLPVars();
    xk_ = std::vector<double>(nx_);
    uk_ = std::vector<double>(nu_);
};

void AdaptiveCorridorHelper::updateNLPVars(){
    nx_ = nlp_->getNx();
    nu_ = nlp_->getNu();
    N_ = nlp_->getN();
    next_index_ = nlp_->getNextInd();
    time_from_index_ = nlp_->getTimeFromIndex();
    final_ind_ = nlp_->getFinalInd();
    free_time_ = nlp_->getFreeTime();

    all_ = std::vector<int>(N_, 0);
    for (int i = 0; i < N_; i++){
        all_[i] = i + (i >= final_ind_);
    }
};

void AdaptiveCorridorHelper::getCollidingPoints(double pos_x, double pos_y, 
                                            double r, std::vector<double>& sol,
                                            std::vector<bool>& constraint_added,
                                            std::vector<int>& to_add,
                                            std::vector<int>& to_remove){
    double margin = 4.0;
    to_add = {}; to_remove = {};

    for (int k = 0; k < N_+1; k++){
        if (k != final_ind_){
            get_x(sol, k);
            if (std::pow(xk_[0] - pos_x, 2) + 
                    std::pow(xk_[1] - pos_y, 2) <= 
                        std::pow(r+margin, 2)){
                if (!constraint_added[k]){
                    to_add.push_back(k);
                    constraint_added[k] = true;
                }
            } else {
                if (constraint_added[k]){
                    to_remove.push_back(k);
                    constraint_added[k] = false;
                }
            }
        }
    }
};

std::vector<int> AdaptiveCorridorHelper::getLeavingPoints(
                                            std::vector<double>& corridor, 
                                            std::vector<double>& sol){
    std::vector<int> k_vals;
    double margin = 2.0;

    for (int k = 0; k < N_+1; k++){
        if (k != final_ind_){
            get_x(sol, k);
            if (xk_[0] >= corridor[1] - margin){
                k_vals.push_back(k);
            }
        }
    }
    return k_vals;
};

std::vector<int> AdaptiveCorridorHelper::getLeavingPoints(
                            std::vector<double>& corridor,
                            std::vector<std::vector<double>>& formatted_xx){
    std::vector<int> k_vals;
    double margin = 2.0;

    for (int k = 0; k < formatted_xx[0].size()-1; k++){
        if (formatted_xx[0][k] >= corridor[1] - margin){
            k_vals.push_back(k);
        }
    }
    return k_vals;
};

bool AdaptiveCorridorHelper::processObstacleConstraints(
                                                    std::vector<double>& x0,
                                                    std::vector<double>& sol, 
                                                    double pos_x, double pos_y,
                                                    double r, int c_idx, 
                                                    int instance_idx){
    if (std::pow(x0[0]-pos_x, 2) + std::pow(x0[1] - pos_y, 2) <= 
            std::pow(r + view_radius_, 2)){
        std::vector<int> to_add = {};
        std::vector<int> to_remove = {};
        getCollidingPoints(pos_x, pos_y, r, sol, 
                           obs_constraint_added_[instance_idx], to_add, 
                           to_remove);
        nlp_->removeExtraConstraint(to_remove, c_idx, instance_idx);
        nlp_->addExtraConstraint(to_add, {c_idx}, {instance_idx}, 
                                 {{pos_x, pos_y, r}});

        return true;
    }
    return false;
};

bool AdaptiveCorridorHelper::processObstacleConstraintsCoarse(
                                                    std::vector<double>& x0,
                                                    double pos_x, double pos_y,
                                                    double r, int c_idx, 
                                                    int instance_idx,
                                                    bool seen){
    if (std::pow(x0[0]-pos_x, 2) + std::pow(x0[1] - pos_y, 2) <= 
            std::pow(r + view_radius_, 2)){
        if (!seen){
            nlp_->addExtraConstraint(all_, {c_idx}, {instance_idx}, 
                                    {{pos_x, pos_y, r}});
        }
        return true;
    }
    nlp_->removeExtraConstraint(all_, c_idx, instance_idx);
    return false;
};

bool AdaptiveCorridorHelper::processLowSpeedZone(std::vector<double>& x0,
                                                    std::vector<double>& sol, 
                                                    double pos_x, double pos_y,
                                                    double r, int c_idx, 
                                                    int instance_idx){
    if (std::pow(x0[0]-pos_x, 2) + std::pow(x0[1] - pos_y, 2) <= 
            std::pow(view_radius_, 2)){
        std::vector<int> to_add = {};
        std::vector<int> to_remove = {};
        getCollidingPoints(pos_x, pos_y, r, sol,
                           low_constraint_added_[instance_idx], to_add, 
                           to_remove);
        nlp_->removeExtraConstraint(to_remove, c_idx, instance_idx);
        nlp_->addExtraConstraint(to_add, {c_idx}, {instance_idx}, 
                                 {{low_speed_min_, r, pos_x, pos_y}});
        return true;
    }
    return false;
};

bool AdaptiveCorridorHelper::processLowSpeedZoneCoarse(std::vector<double>& x0,
                                                    double pos_x, double pos_y,
                                                    double r, int c_idx, 
                                                    int instance_idx, 
                                                    bool seen){
    if (std::pow(x0[0]-pos_x, 2) + std::pow(x0[1] - pos_y, 2) <= 
            std::pow(view_radius_, 2)){
        if (!seen){
            nlp_->addExtraConstraint(all_, {c_idx}, {instance_idx}, 
                                    {{low_speed_min_, r, pos_x, pos_y}});
        }
        return true;
    }
    nlp_->removeExtraConstraint(all_, c_idx, instance_idx);
    return false;
};

std::vector<double> AdaptiveCorridorHelper::shiftInitialGuess(
                                                    std::vector<double>& sol){
    std::vector<double> res(sol.size(), 0.0);

    int k = 0;
    int next_k = next_index_[k].value();
    while (next_k != final_ind_){
        get_x(sol, next_k);
        std::copy(xk_.begin(), xk_.end(), res.begin() + 
            (nx_+nu_)*k - nu_*(k > final_ind_));
        get_u(sol, next_k);
        std::copy(uk_.begin(), uk_.end(), res.begin() + 
            (nx_+nu_)*k + nx_ - nu_*(k > final_ind_));
        k = next_k;
        next_k = next_index_[next_k].value();
    }

    get_x(sol, final_ind_);
    std::copy(xk_.begin(), xk_.end(), res.begin() + 
        (nx_+nu_)*k - nu_*(k > final_ind_));
    get_u(sol, k);
    std::copy(uk_.begin(), uk_.end(), res.begin() + 
        (nx_+nu_)*k + nx_ - nu_*(k > final_ind_));

    std::copy(xk_.begin(), xk_.end(), res.begin() + 
        (nx_+nu_)*final_ind_ - nu_*(final_ind_ > final_ind_));
    
    return res;
};

void AdaptiveCorridorHelper::appendInitGuess(std::vector<double>& sol, 
                                             std::vector<double>& init_guess,
                                             int nb_steps){
    get_x(sol, final_ind_);
    int k = 0;
    while (next_index_[k] != final_ind_){
        k = next_index_[k].value();
    }
    get_u(sol, k);
    for (int i = 0; i < nb_steps-1; i++){
        if (i == 0){
            for (int j = 0; j < nu_; j++){ init_guess.push_back(uk_[j]);}
        }
        for (int j = 0; j < nx_; j++){ init_guess.push_back(xk_[j]);}
        if (i < nb_steps-2){
            for (int j = 0; j < nu_; j++){ init_guess.push_back(uk_[j]);}
        }
    }
};

std::vector<std::vector<double>> AdaptiveCorridorHelper::performAdaptiveLoop(
        AdaptiveNLP& adaptiveNLP, Plotter& plotter, std::vector<double> x0, 
        std::vector<double>& corridor_1, std::vector<double>& corridor_2, 
        int nx, int nu, int nb_steps, int nb_iterations, double dt, 
        int nb_intervals, std::vector<double>& obstacle_x, 
        std::vector<double>& obstacle_y, std::vector<double>& obstacle_r, 
        std::vector<double>& low_speed_pos_x, 
        std::vector<double>& low_speed_pos_y, 
        std::vector<double>& low_speed_r){
    
    // remove ipopt printing
    adaptiveNLP.setIpoptPrintLevel(0);

    ////////////////////////////////////////
    // define containers to store results //
    ////////////////////////////////////////

    // for every MPC iteration, store which obstacles/people are visible
    std::vector<std::vector<bool>> visible_constraints(nb_iterations,
        std::vector<bool>(obstacle_x.size() + low_speed_pos_x.size(), false));

    // solution times
    std::vector<double> timings(nb_iterations);

    // total times
    std::vector<double> timings_total(nb_iterations);

    // initial states (traced path)
    std::vector<std::vector<double>> x0s (nb_iterations, 
                                          std::vector<double>(x0.size()));

    // number of constraints present in the NLP
    std::vector<int> nb_constraints(nb_iterations);
    std::vector<std::vector<int>> nb_constraints_specific(nb_iterations, 
        std::vector<int>(3+3)); // boundary - fixed - discretization - extra
    
    // computed planned trajectories
    std::vector<std::vector<std::vector<double>>> solutions(nb_iterations,
        std::vector<std::vector<double>>(nx+nu));

    // number of solver iterations
    std::vector<int> iter_counts(nb_iterations);

    // variable to store IPOPT solution time
    double computation_time;

    // variable to store planned state trajectories in a formatted way
    std::vector<std::vector<double>> formatted_xx;

    // variable to store planned control trajectories in a formatted way
    std::vector<std::vector<double>> formatted_uu;

    // get bookkeeping info about chronological sequence of time-steps
    std::vector<std::optional<int>> next_ind = 
        adaptiveNLP.getNextInd();

    // get time-grid
    std::vector<std::optional<double>> time_from_index = 
        adaptiveNLP.getTimeFromIndex();

    // update plotter with NLP information (so it can properly format 
    // solutions)
    plotter.update(adaptiveNLP.getN(), adaptiveNLP.getFinalInd(), 
                    next_ind, time_from_index);

    obs_constraint_added_ = std::vector<std::vector<bool>>(adaptiveNLP.getN(),
        std::vector<bool>(obstacle_x.size(), false));
    low_constraint_added_ = std::vector<std::vector<bool>>(adaptiveNLP.getN(),
        std::vector<bool>(low_speed_pos_x.size(), false));

    //////////////////////////////////
    // Solve NLP for the first time //
    //////////////////////////////////
    adaptiveNLP.solveNlp({{"p_g0", x0}, {"p_gT", corridor_1}}, 
                         computation_time);

    // store solution
    x0s[0] = x0;
    std::vector<double> sol = adaptiveNLP.getSolution();
    plotter.formatSolution(sol, formatted_xx, formatted_uu);
    for (int i = 0; i < nx; i++){
        solutions[0][i] = formatted_xx[i];
    }
    for (int i = 0; i < nu; i++){
        solutions[0][nx+i] = formatted_uu[i];
    }
    timings[0] = computation_time;
    timings_total[0] = computation_time;
    nb_constraints[0] = adaptiveNLP.getNumberOfConstraints();
    nb_constraints_specific[0][0] = 
        adaptiveNLP.getNumberOfBoundaryConstraints();
    nb_constraints_specific[0][1] = nb_constraints_specific[0][0] +
        adaptiveNLP.getNumberOfFixedConstraints();
    nb_constraints_specific[0][2] = nb_constraints_specific[0][1] +
        adaptiveNLP.getNumberOfDiscretizationConstraints();
    std::vector<int> nb_extras = adaptiveNLP.getNumberOfExtraConstraints();
    for (int i = 0; i < nb_extras.size(); i++){
        nb_constraints_specific[0][3+i] = nb_constraints_specific[0][2+i] +
            nb_extras[i];
    }
    iter_counts[0] = adaptiveNLP.getIterCount();

    ////////////////////////////////////////////////////////////////////
    // allocate some vectors to avoid allocations during the MPC loop //
    ////////////////////////////////////////////////////////////////////
    
    // vector containing right-hand-side of equation x_dot = f(x, u)
    std::vector<double> rhs(nx);

    // allocated time-steps
    std::vector<int> k_vals;

    // initial guess
    std::vector<double> init_guess;

    // corridor to which the terminal state belongs
    std::vector<double> final_point_corridor = corridor_1;

    // bool indicating visibility
    bool v;

    // number of time steps with which to append to horizon in the case of 
    // starting with a shorter horizon
    int nb_steps_to_add = 2;
    std::vector<int> extraSteps(nb_steps_to_add-1);

    ////////////////////////////////
    // MPC loop using AdaptiveNLP //
    ////////////////////////////////

    // bool indicating whether to add no-collision constraints to only a 
    // subset of time-steps
    const bool USE_FINE_APPROACH = true;

    // for every MPC iteration
    for (int iteration = 1; iteration < nb_iterations; iteration++){
        updateNLPVars();
        std::cout<<"Iteration "<<iteration<<endl;
        auto start = high_resolution_clock::now();

        // update initial position
        for (int i = 0; i < nx; i++){ x0[i] = sol[(nx+nu) + i];}
        x0s[iteration] = x0;

        // remove points from corridor 1 and add them to corridor 2
        k_vals = getLeavingPoints(corridor_1, sol);
        if (k_vals.size() > 0){ final_point_corridor = corridor_2;};
        adaptiveNLP.removeExtraConstraint(k_vals, 0);
        adaptiveNLP.addExtraConstraint(k_vals, {0}, {corridor_2});

        // process all obstacles
        for (int i = 0; i < obstacle_x.size(); i++){
            if (USE_FINE_APPROACH){
                v = processObstacleConstraints(x0, sol, obstacle_x[i], 
                                               obstacle_y[i], obstacle_r[i], 
                                               1, i);
            } else {
                v = processObstacleConstraintsCoarse(x0, obstacle_x[i], 
                                                     obstacle_y[i], 
                                                     obstacle_r[i], 1, i, 
                                                     visible_constraints
                                                        [iteration-1][i]);
            }
            visible_constraints[iteration][i] = v;
        }

        // process low-speed zones
        for (int i = 0; i < low_speed_pos_x.size(); i++){
            if (USE_FINE_APPROACH){
                v = processLowSpeedZone(x0, sol, low_speed_pos_x[i], 
                                        low_speed_pos_y[i], low_speed_r[i], 
                                        2, i);
            } else {
                v = processLowSpeedZoneCoarse(x0, low_speed_pos_x[i], 
                                              low_speed_pos_y[i], 
                                              low_speed_r[i], 2, i, 
                                              visible_constraints[iteration-1]
                                                [obstacle_x.size() + i]);
            }
            visible_constraints[iteration][obstacle_x.size() + i] = v;
        }

        // remove points from corridor 2
        k_vals = getLeavingPoints(corridor_2, sol);
        if (k_vals.size() > 0){ 
            final_point_corridor = {-1000.0, 1000.0, -1000.0, 1000.0};
        };
        adaptiveNLP.removeExtraConstraint(k_vals, 0);

        // get new initial guess
        init_guess = shiftInitialGuess(sol);

        // clear structure
        adaptiveNLP.clearStructuralZeros();

        // potenitally append the horizon
        if (iteration < 10 && adaptiveNLP.getN() < nb_intervals*(nb_steps-1)){
            // update NLP
            adaptiveNLP.appendTimeSteps({nb_steps_to_add});

            // update Corridor constraints
            for (int i = 0; i < nb_steps_to_add-1; i++){
                extraSteps[i] = adaptiveNLP.getN() - (nb_steps_to_add-1) + 1 + i;
            }

            adaptiveNLP.addExtraConstraint(extraSteps, {0}, {corridor_1});

            // update initial guess
            appendInitGuess(sol, init_guess, nb_steps_to_add);
        }
        
        // Solve the NLP
        cout<<"EXIT: "<<adaptiveNLP.solveNlp({{"p_g0", x0}, 
                                              {"p_gT", final_point_corridor}}, 
                                             computation_time, init_guess)<<
                                             endl;

        auto stop = high_resolution_clock::now();

        // store logging information
        x0s[iteration] = x0;
        sol = adaptiveNLP.getSolution();
        timings[iteration] = computation_time;
        timings_total[iteration] = 0*computation_time +
            double(duration_cast<nanoseconds>(stop-start).count())*1.0e-6;
        nb_constraints[iteration] = adaptiveNLP.getNumberOfConstraints();
        nb_constraints_specific[iteration][0] = 
            adaptiveNLP.getNumberOfBoundaryConstraints();
        nb_constraints_specific[iteration][1] = 
            nb_constraints_specific[iteration][0] +
            adaptiveNLP.getNumberOfFixedConstraints();
        nb_constraints_specific[iteration][2] = 
            nb_constraints_specific[iteration][1] +
            adaptiveNLP.getNumberOfDiscretizationConstraints();
        std::vector<int> nb_extras = adaptiveNLP.getNumberOfExtraConstraints();
        for (int i = 0; i < nb_extras.size(); i++){
            nb_constraints_specific[iteration][3+i] = 
                nb_constraints_specific[iteration][2+i] +
                nb_extras[i];
        }
        iter_counts[iteration] = adaptiveNLP.getIterCount();

        std::vector<std::optional<int>> next_ind = 
            adaptiveNLP.getNextInd();
        std::vector<std::optional<double>> time_from_index = 
            adaptiveNLP.getTimeFromIndex();
        plotter.update(adaptiveNLP.getN(), adaptiveNLP.getFinalInd(), 
                        next_ind, time_from_index);
        plotter.formatSolution(sol, formatted_xx, formatted_uu);
            for (int i = 0; i < nx; i++){
            solutions[iteration][i] = formatted_xx[i];
        }
        for (int i = 0; i < nu; i++){
            solutions[iteration][nx+i] = formatted_uu[i];
        }
    }

    // store all logging information with the plotter object
    plotter.addAdaptiveResults(solutions, x0s, visible_constraints, 
                            nb_constraints, nb_constraints_specific, 
                            timings, timings_total);
    plotter.addAdaptiveIterCounts(iter_counts);

    return {timings, timings_total};
};

std::vector<std::vector<double>> AdaptiveCorridorHelper::performCasadiLoop(
            BuildingBlocks& blocks, Plotter& plotter, 
            std::vector<double> x0, std::vector<double>& corridor_1, 
            std::vector<double>& corridor_2, int nx, int nu, int nb_steps, 
            int nb_iterations, double dt, double T, int nb_intervals, 
            std::vector<double>& obstacle_x, std::vector<double>& obstacle_y, 
            std::vector<double>& obstacle_r, 
            std::vector<double>& low_speed_pos_x, 
            std::vector<double>& low_speed_pos_y, 
            std::vector<double>& low_speed_r){
    
    int N = nb_intervals*(nb_steps-1);
    
    ////////////////////////////////////////
    // define containers to store results //
    ////////////////////////////////////////

    // storage for formatted state and control trajectories
    std::vector<std::vector<double>> formatted_xx(nx, 
        std::vector<double>(N+1));
    std::vector<std::vector<double>> formatted_uu(nu, 
        std::vector<double>(N));

    // for every MPC iteration, store which obstacles/people are visible
    std::vector<std::vector<bool>> visible_constraints(nb_iterations,
        std::vector<bool>(obstacle_x.size() + low_speed_pos_x.size(), false));

    // solution times
    std::vector<double> timings(nb_iterations);

    // total times
    std::vector<double> timings_total(nb_iterations);

    // initial states (traced path)
    std::vector<std::vector<double>> x0s (nb_iterations, 
                                          std::vector<double>(x0.size()));

    // number of constraints present in the NLP
    std::vector<int> nb_constraints(nb_iterations);
    std::vector<std::vector<std::vector<double>>> solutions(nb_iterations,
        std::vector<std::vector<double>>(nx+nu)); 

    // number of solver iterations
    std::vector<int> iter_counts(nb_iterations);
    
    ////////////////////////////////////
    // construct CasADi Opti instance //
    ////////////////////////////////////
    Opti opti = Opti();
    // states
    MX xx = opti.variable(nx, N+1);
    // controls
    MX uu = opti.variable(nu, N);
    // initial state constraint parameter
    MX p_0 = opti.parameter(nx, 1);
    // terminal state constraint parameter
    MX p_T = opti.parameter(4, 1);
    // corridor constraint parameter
    MX p_corridors = opti.parameter(4, N);
    // no-collision constraint parameter
    MX p_obs = opti.parameter(3, obstacle_x.size());
    // safety constraint parameter
    MX p_low = opti.parameter(4, low_speed_pos_x.size());

    // containers to store initial guesses
    DM xx_init = DM(xx.size());
    DM uu_init = DM(uu.size());

    // define objective
    MX obj = blocks.eval_Phi_0({xx(Slice(),0), 0})[0] + 
             blocks.eval_Phi_f({xx(Slice(),N), 0})[0];
    for (int k = 0; k < N; k++){
        obj += blocks.eval_phi({xx(Slice(), k), uu(Slice(), k), dt, 0})[0];
    }
    opti.minimize(obj);

    // define constraints
    // initial and terminal state constraints
    opti.subject_to((blocks.get_g0_lb() <= 
                    blocks.eval_g0({xx(Slice(),0), p_0})[0]) <= 
                    blocks.get_g0_ub());
    opti.subject_to((blocks.get_gT_lb() <= 
                    blocks.eval_gT({xx(Slice(),N), p_T})[0]) <= 
                    blocks.get_gT_ub());

    // fixed and extra constraint
    for (int k = 0; k < N; k++){
        opti.subject_to((blocks.get_g_fixed_lb() <= 
            blocks.eval_g_fixed({xx(Slice(), k), uu(Slice(), k), 0})[0]) <= 
            blocks.get_g_fixed_ub());
        opti.subject_to((blocks.get_g_extra_lb(0) <= 
            blocks.eval_g_extra({xx(Slice(), k), uu(Slice(), k), 
                                p_corridors(Slice(), k)}, 0)[0]) <= 
            blocks.get_g_extra_ub(0));
        
        for (int i = 0; i < obstacle_x.size(); i++){
            opti.subject_to((blocks.get_g_extra_lb(1) <= 
                blocks.eval_g_extra({xx(Slice(), k), uu(Slice(), k), 
                                    p_obs(Slice(), i)}, 1)[0]) <= 
                blocks.get_g_extra_ub(1));
        }
        for (int i = 0; i < low_speed_pos_x.size(); i++){
            opti.subject_to((blocks.get_g_extra_lb(2) <= 
                blocks.eval_g_extra({xx(Slice(), k), uu(Slice(), k), 
                                    p_low(Slice(), i)}, 2)[0]) <= 
                blocks.get_g_extra_ub(2));
        }
    }

    // dynamics
    assert (nb_steps == 2 || nb_steps == 3);
    for (int j = 0; j < nb_intervals; j++){
        if (nb_steps == 2){
            MX local_xx = horzcat(xx(Slice(), j), xx(Slice(), j+1));
            opti.subject_to(
                blocks.eval_g_disc
                    ({local_xx, uu(Slice(), j), dt, T}, nb_steps)[0] == 0);
        } else {
            MX local_xx = horzcat(xx(Slice(), 2*j), xx(Slice(), 2*j+1),
                                  xx(Slice(), 2*j+2));
            MX local_uu = horzcat(uu(Slice(), 2*j), uu(Slice(), 2*j+1));
            opti.subject_to(
                blocks.eval_g_disc
                    ({local_xx, local_uu, 
                        std::vector<double>{dt, dt}, T}, nb_steps)[0] == 0);
        }
    }

    // set initial guess
    for (int i = 0; i < N+1; i++){ 
        xx_init(Slice(), i) = blocks.eval_x_init({DM(dt*i)})[0];
    }
    for (int i = 0; i < N; i++){ 
        uu_init(Slice(), i) = blocks.eval_u_init({DM(dt*i)})[0];
    }

    // set parameters
    opti.solver("ipopt", {}, {{"print_level", 0}});

    // set parameter values
    DM p_0_DM = DM(x0);
    DM p_T_DM = DM(corridor_1);
    DM p_corridors_DM = DM(p_corridors.size1(), p_corridors.size2());
    DM p_obs_DM = DM(p_obs.size());
    DM p_low_DM = DM(p_low.size());

    for (int k = 0; k < N; k++){
        p_corridors_DM(Slice(), k) = DM(corridor_1);
    }
    for (int i = 0; i < obstacle_x.size(); i++){
        // define dummy obstacle far away
        p_obs_DM(Slice(), i) = DM({1000.0+1000.0*i, 1000.0+1000.0*i, 0.1});
    }
    for (int i = 0; i < low_speed_pos_x.size(); i++){
        // define dummy people far away
        p_low_DM(Slice(), i) = DM({low_speed_min_, low_speed_r[i], 
                       1000.0+1000.0*i, 1000.0+1000.0*i});
    }
    
    // create opti function object (more efficient)
    Function opti_function = opti.to_function("opti_function", 
        {p_0, p_T, p_corridors, p_obs, p_low, opti.x()}, {xx, uu}, 
        {{"record_time", true}});

    // define opti function input vector
    std::vector<DM> args = {p_0_DM, p_T_DM, p_corridors_DM, p_obs_DM, 
        p_low_DM, 
        vertcat(reshape(xx_init, nx*(N+1), 1), reshape(uu_init, nu*N, 1))};

    //////////////////////////////////
    // Solve NLP for the first time //
    //////////////////////////////////
    std::vector<DM> output = opti_function(args);

    // store solution
    DM xx_DM = output[0];
    DM uu_DM = output[1];

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
    solutions[0][0] = formatted_xx[0];
    solutions[0][1] = formatted_xx[1];
    solutions[0][2] = formatted_xx[2];
    solutions[0][3] = formatted_xx[3];
    solutions[0][4] = formatted_uu[0];
    solutions[0][5] = formatted_uu[1];

    timings[0] = 1000*double(opti_function.stats()["t_wall_total"]);
    timings_total[0] = timings[0];
    x0s[0] = x0;
    nb_constraints[0] = 4 + 4 + N*(4 + 4 + 4 + obstacle_x.size()*1 + 
                                                    low_speed_pos_x.size()*1);
    iter_counts[0] = 0;
    
    ////////////////////////////////////////////////////////////////////
    // allocate some vectors to avoid allocations during the MPC loop //
    ////////////////////////////////////////////////////////////////////

    // vector containing right-hand-side of equation x_dot = f(x, u)
    std::vector<double> rhs(nx);
    
    // allocated time-steps
    std::vector<int> k_vals;
    
    // initial guess
    std::vector<double> init_guess;
    
    // bool indicating visibility
    bool v;
    
    // store to which corridor every point belongs
    std::vector<int> corridor_id(N, 1);

    // define dummy corridor
    std::vector<double> final_corridor = {-1000.0, 1000.0, -1000.0, 1000.0};

    /////////////////////////////////////
    // MPC loop for case CasADi Opti 1 //
    /////////////////////////////////////
    for (int iteration = 1; iteration < nb_iterations; iteration++){
        std::cout<<"Iteration "<<iteration<<endl;
        auto start = high_resolution_clock::now();

        // update initial position
        for (int i = 0; i < nx; i++){ x0[i] = formatted_xx[i][1];}
        x0s[iteration] = x0;
        p_0_DM = DM(x0);

        // remove points from corridor 1 and add them to corridor 2
        k_vals = getLeavingPoints(corridor_1, formatted_xx);
        for (int k : k_vals){
            if (corridor_id[k] == 1){
                corridor_id[k] = 2;
                p_corridors_DM(Slice(), k) = DM(corridor_2);
            }
        }

        // process obstacles
        for (int i = 0; i < obstacle_x.size(); i++){
            v = false;
            if (std::pow(obstacle_x[i]-x0[0], 2) + 
                    std::pow(obstacle_y[i]-x0[1], 2) <= 
                    std::pow(obstacle_r[i]+view_radius_, 2)){
                p_obs_DM(Slice(), i) = DM({obstacle_x[i], obstacle_y[i], 
                                          obstacle_r[i]});
                v = true;
            }
            visible_constraints[iteration][i] = v;
        }

        // process low-speed zones
        for (int i = 0; i < low_speed_pos_x.size(); i++){
            v = false;
            if (std::pow(low_speed_pos_x[i]-x0[0], 2) + 
                    std::pow(low_speed_pos_y[i]-x0[1], 2) <= 
                    std::pow(view_radius_, 2)){
                p_low_DM(Slice(), i) = DM({low_speed_min_, low_speed_r[i],
                               low_speed_pos_x[i], low_speed_pos_y[i]});
                v = true;
            }
            visible_constraints[iteration][obstacle_x.size()+i] = v;
        }

        // remove points from corridor 2
        k_vals = getLeavingPoints(corridor_2, formatted_xx);
        for (int k : k_vals){
            if (corridor_id[k] == 2){
                p_corridors_DM(Slice(), k) = DM(final_corridor);
                corridor_id[k] = 3;
            }
        } 

        // create new initial guess
        xx_init(Slice(), Slice(0, N-1)) = xx_DM(Slice(), Slice(1, N));
        xx_init(Slice(), N-1) = xx_DM(Slice(), N);
        xx(Slice(), N) = xx_DM(Slice(), N);
        uu(Slice(), Slice(0, N-1)) = uu_DM(Slice(), Slice(1, N));
        uu(Slice(), N-1) = uu_DM(Slice(), N-1);

        if (corridor_id[corridor_id.size()-1] == 1){
            p_T_DM = DM(corridor_1);
        } else if (corridor_id[corridor_id.size()-1] == 2){
            p_T_DM = DM(corridor_2);
        } else {
            p_T_DM = DM(final_corridor);
        }

        // Solve NLP
        args = {p_0_DM, p_T_DM, p_corridors_DM, p_obs_DM, p_low_DM, 
            vertcat(reshape(xx_init, nx*(N+1), 1), reshape(uu_init, nu*N, 1))};
        output = opti_function(args);

        // Store output
        xx_DM = output[0];
        uu_DM = output[1];
        auto stop = high_resolution_clock::now();
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
        solutions[iteration][0] = formatted_xx[0];
        solutions[iteration][1] = formatted_xx[1];
        solutions[iteration][2] = formatted_xx[2];
        solutions[iteration][3] = formatted_xx[3];
        solutions[iteration][4] = formatted_uu[0];
        solutions[iteration][5] = formatted_uu[1];

        timings[iteration] = 1000*double(opti_function.stats()["t_wall_total"]);
        timings_total[iteration] = 0*timings[iteration] + 
            double(duration_cast<nanoseconds>(stop-start).count())*1.0e-6;
        x0s[iteration] = x0;
        nb_constraints[iteration] = 4 + 4 + 
            N*(4 + 4 + 4 + obstacle_x.size()*1 + low_speed_pos_x.size()*1);
        iter_counts[iteration] = 0;
    }

    // store all logging information with the plotter object
    plotter.addRefResults(solutions, x0s, visible_constraints, 
                            nb_constraints, {{}}, timings, timings_total);
    plotter.addRefIterCounts(iter_counts);

    return {timings, timings_total};
};

std::vector<std::vector<double>> AdaptiveCorridorHelper::performCasadiLoopReformulation(
        BuildingBlocks& blocks, Plotter& plotter, 
        std::vector<double> x0, std::vector<double>& corridor_1, 
        std::vector<double>& corridor_2, int nx, int nu, int nb_steps, 
        int nb_iterations, double dt, double T, int nb_intervals, 
        std::vector<double>& obstacle_x, std::vector<double>& obstacle_y, 
        std::vector<double>& obstacle_r, 
        std::vector<double>& low_speed_pos_x, 
        std::vector<double>& low_speed_pos_y, 
        std::vector<double>& low_speed_r){
    int N = nb_intervals*(nb_steps-1);
    
    ////////////////////////////////////////
    // define containers to store results //
    ////////////////////////////////////////

    // storage for formatted state and control trajectories
    std::vector<std::vector<double>> formatted_xx(nx, 
        std::vector<double>(N+1));
    std::vector<std::vector<double>> formatted_uu(nu, 
        std::vector<double>(N));

    // for every MPC iteration, store which obstacles/people are visible
    std::vector<std::vector<bool>> visible_constraints(nb_iterations,
        std::vector<bool>(obstacle_x.size() + low_speed_pos_x.size(), false));

    // solution times
    std::vector<double> timings(nb_iterations);

    // total times
    std::vector<double> timings_total(nb_iterations);\

    // initial states (traced path)
    std::vector<std::vector<double>> x0s (nb_iterations, 
                                          std::vector<double>(x0.size()));

    // number of constraints present in the NLP
    std::vector<int> nb_constraints(nb_iterations);
    std::vector<std::vector<std::vector<double>>> solutions(nb_iterations,
        std::vector<std::vector<double>>(nx+nu)); 
    
    // number of solver iterations
    std::vector<int> iter_counts(nb_iterations);
    
    ////////////////////////////////////
    // construct CasADi Opti instance //
    ////////////////////////////////////
    double corridor_margin = 2;
    double obstacle_margin = 4;
    double low_speed_margin = 4;
    int constraint_counter;
    std::vector<double> p_obs(3);
    std::vector<double> p_low(4);
    std::vector<double> init_guess;
    std::vector<double> final_point_corridor = corridor_1;
    bool v;
    DM xx_DM;
    DM uu_DM;

    /////////////////////////////////////
    // MPC loop for case CasADi Opti 2 //
    /////////////////////////////////////
    for (int iteration = 0; iteration < nb_iterations; iteration++){
        cout<<"Iteration "<<iteration<<endl;
        auto start = high_resolution_clock::now();
        constraint_counter = 0;
        x0s[iteration] = x0;

        // create opti instance
        Opti opti = Opti();

        // states
        MX xx = opti.variable(nx, N+1);
        // controls
        MX uu = opti.variable(nu, N);

        // objective
        MX obj = blocks.eval_Phi_0({xx(Slice(),0), 0})[0] + 
                blocks.eval_Phi_f({xx(Slice(),N), 0})[0];
        for (int k = 0; k < N; k++){
            obj += blocks.eval_phi({xx(Slice(), k), uu(Slice(), k), dt, 0})[0];
        }
        opti.minimize(obj);

        // constraints
        // initial state constraint
        opti.subject_to((blocks.get_g0_lb() <= 
                        blocks.eval_g0({xx(Slice(),0), x0})[0]) <= 
                        blocks.get_g0_ub());
        constraint_counter += blocks.get_nb_g0();
        for (int k = 0; k < N; k++){
            // fixed constraints
            opti.subject_to((blocks.get_g_fixed_lb() <= 
                blocks.eval_g_fixed({xx(Slice(), k), uu(Slice(), k), 0})[0]) <= 
                blocks.get_g_fixed_ub());
            constraint_counter += blocks.get_nb_g_fixed();

            // corridor constraint
            if (iteration > 0){
                if (formatted_xx[0][k] < corridor_2[1] - corridor_margin){
                    if (formatted_xx[0][k] < corridor_1[1] - corridor_margin){
                        opti.subject_to((blocks.get_g_extra_lb(0) <= 
                            blocks.eval_g_extra({xx(Slice(), k), uu(Slice(), k), 
                                        corridor_1}, 0)[0]) <= 
                            blocks.get_g_extra_ub(0));
                        final_point_corridor = corridor_1;
                        constraint_counter += blocks.get_nb_g_extra(0);
                    } else {
                        opti.subject_to((blocks.get_g_extra_lb(0) <= 
                            blocks.eval_g_extra({xx(Slice(), k), uu(Slice(), k), 
                                        corridor_2}, 0)[0]) <= 
                            blocks.get_g_extra_ub(0));
                        constraint_counter += blocks.get_nb_g_extra(0);
                        final_point_corridor = corridor_2;
                    }
                } else {
                    final_point_corridor = {-1000, 1000, -1000, 1000};
                }
            } else {
                opti.subject_to((blocks.get_g_extra_lb(0) <= 
                    blocks.eval_g_extra({xx(Slice(), k), uu(Slice(), k), 
                                        corridor_1}, 0)[0]) <= 
                    blocks.get_g_extra_ub(0));
                constraint_counter += blocks.get_nb_g_extra(0);
            }
            
            // obstacle constraint
            if (iteration > 0){
                for (int i = 0; i < obstacle_x.size(); i++){
                    v = false;
                    if (std::pow(x0[0] - obstacle_x[i], 2) + 
                        std::pow(x0[1] - obstacle_y[i], 2) <= 
                            std::pow(obstacle_r[i] + view_radius_, 2) &&
                        std::pow(formatted_xx[0][k] - obstacle_x[i], 2) + 
                        std::pow(formatted_xx[1][k] - obstacle_y[i], 2) <= 
                            std::pow(obstacle_r[i] + obstacle_margin, 2)){
                        p_obs = {obstacle_x[i], obstacle_y[i], obstacle_r[i]};
                        opti.subject_to((blocks.get_g_extra_lb(1) <= 
                            blocks.eval_g_extra({xx(Slice(), k), uu(Slice(), k), 
                                p_obs}, 1)[0]) <= 
                            blocks.get_g_extra_ub(1));
                        constraint_counter += blocks.get_nb_g_extra(1);
                        v = true;
                    };
                    visible_constraints[iteration][i] = v;
                }
            }

            // low speed constraint
            if (iteration > 0){
                for (int i = 0; i < low_speed_pos_x.size(); i++){
                    v = false;
                    if (std::pow(x0[0] - low_speed_pos_x[i], 2) + 
                        std::pow(x0[1] - low_speed_pos_y[i], 2) <= 
                            std::pow(view_radius_, 2) &&
                        std::pow(formatted_xx[0][k] - low_speed_pos_x[i], 2) + 
                        std::pow(formatted_xx[1][k] - low_speed_pos_y[i], 2) <= 
                            std::pow(low_speed_r[i] + low_speed_margin, 2)){
                        p_low = {low_speed_min_, low_speed_r[i], 
                                 low_speed_pos_x[i], low_speed_pos_y[i]};
                        opti.subject_to((blocks.get_g_extra_lb(2) <= 
                            blocks.eval_g_extra({xx(Slice(), k), uu(Slice(), k), 
                                p_low}, 2)[0]) <= 
                            blocks.get_g_extra_ub(2));
                        constraint_counter += blocks.get_nb_g_extra(2);
                        v = true;
                    };
                    visible_constraints[iteration][obstacle_x.size() + i] = v;
                }
            }
        }

        // terminal state constraint
        opti.subject_to((blocks.get_gT_lb() <= 
                        blocks.eval_gT({xx(Slice(),N), final_point_corridor})[0]) 
                        <= blocks.get_gT_ub());
        constraint_counter += blocks.get_nb_gT();

        assert (nb_steps == 2 || nb_steps == 3);
        for (int j = 0; j < nb_intervals; j++){
            if (nb_steps == 2){
                MX local_xx = horzcat(xx(Slice(), j), xx(Slice(), j+1));
                opti.subject_to(
                    blocks.eval_g_disc
                        ({local_xx, uu(Slice(), j), dt, T}, nb_steps)[0] == 0);
                constraint_counter += blocks.get_nb_g_disc(nb_steps);
            } else {
                MX local_xx = horzcat(xx(Slice(), 2*j), xx(Slice(), 2*j+1),
                                    xx(Slice(), 2*j+2));
                MX local_uu = horzcat(uu(Slice(), 2*j), uu(Slice(), 2*j+1));
                opti.subject_to(
                    blocks.eval_g_disc
                        ({local_xx, local_uu, 
                            std::vector<double>{dt, dt}, T}, nb_steps)[0] == 0);
                constraint_counter += blocks.get_nb_g_disc(nb_steps);
            }
        }

        // initial guess
        if (iteration > 0){
            for (int k = 0; k < N-1; k++){
                opti.set_initial(xx(Slice(), k), xx_DM(Slice(), k+1));
                opti.set_initial(uu(Slice(), k), uu_DM(Slice(), k+1));
            }
            opti.set_initial(xx(Slice(), N-1), xx_DM(Slice(), N));
            opti.set_initial(xx(Slice(), N), xx_DM(Slice(), N));
            opti.set_initial(uu(Slice(), N-1), uu_DM(Slice(), N-1));
        } else {
            for (int i = 0; i < N+1; i++){ 
                opti.set_initial(xx(Slice(), i), 
                                 blocks.eval_x_init({DM(dt*i)})[0]);
            }
            for (int i = 0; i < N; i++){ 
                opti.set_initial(uu(Slice(), i), 
                                 blocks.eval_u_init({DM(dt*i)})[0]);
            }
        }

        // set parameters
        opti.solver("ipopt", {}, {{"print_level", 0}});
        
        // Solve NLP
        OptiSol sol_casadi = opti.solve();

        // Store output
        xx_DM = sol_casadi.value(xx);
        uu_DM = sol_casadi.value(uu);
        auto stop = high_resolution_clock::now();
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
        solutions[iteration][0] = formatted_xx[0];
        solutions[iteration][1] = formatted_xx[1];
        solutions[iteration][2] = formatted_xx[2];
        solutions[iteration][3] = formatted_xx[3];
        solutions[iteration][4] = formatted_uu[0];
        solutions[iteration][5] = formatted_uu[1];

        timings[iteration] = 1000.0*double(sol_casadi.stats()["t_wall_total"]);
        timings_total[iteration] = 0*timings[iteration] + 
            double(duration_cast<nanoseconds>(stop-start).count())*1.0e-6;
        for (int i = 0; i < nx; i++){x0[i] = formatted_xx[i][1];};
        nb_constraints[iteration] = constraint_counter;
        iter_counts[iteration] = sol_casadi.stats()["iter_count"];
    }

    // store all logging information with the plotter object
    plotter.addNaiveResults(solutions, x0s, visible_constraints, 
                            nb_constraints, {{}}, timings, timings_total);
    plotter.addNaiveIterCounts(iter_counts);

    return {timings, timings_total};
};

void AdaptiveCorridorHelper::get_x(std::vector<double>& sol, int k){
	int offset = (k > final_ind_)*nu_;
	for (int i = 0; i < nx_; i++){
		xk_[i] = sol[(nx_ + nu_)*k - offset + free_time_ + i];
	}
}

void AdaptiveCorridorHelper::get_u(std::vector<double>& sol, int k){
	int offset = (k > final_ind_)*nu_;
	for (int i = 0; i < nu_; i++){
		uk_[i] = sol[(nx_ + nu_)*k + nx_ - offset + free_time_ + i];
	}
}