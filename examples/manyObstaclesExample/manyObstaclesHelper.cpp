#include "manyObstaclesHelper.hpp"
#include <chrono>

using namespace chrono;

ManyObstaclesHelper::ManyObstaclesHelper(AdaptiveNLP& nlp):nlp_(&nlp){
    updateNLPVars();
    xk_ = std::vector<double>(nx_);
    uk_ = std::vector<double>(nu_);
};

void ManyObstaclesHelper::updateNLPVars(){
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

std::vector<int> ManyObstaclesHelper::getCollidingPoints(double pos_x, 
                                                    double pos_y, double r, 
                                                    std::vector<double>& sol,
                                                    double margin){
    std::vector<int> k_vals;

    for (int k = 0; k < N_+1; k++){
        if (k != final_ind_){
            get_x(sol, k);
            if (std::pow(xk_[0] - pos_x, 2) + 
                    std::pow(xk_[1] - pos_y, 2) <= 
                        std::pow(r+margin, 2)){
                k_vals.push_back(k);
            }
        }
    }
    return k_vals;
};

bool ManyObstaclesHelper::checkFeasible(std::vector<double>& sol,
                                        std::vector<double>& xx,
                                        std::vector<double>& yy,
                                        std::vector<double>& rr){
    std::vector<int> k_vals;
    for (int i = 0; i < xx.size(); i++){
        k_vals = getCollidingPoints(xx[i], yy[i], rr[i], sol, -1.0e-6);
        if (k_vals.size() > 0){ return false;}
    }
    return true;
};

bool ManyObstaclesHelper::processObstacleConstraints(std::vector<double>& sol, 
                                                     std::vector<double> xx, 
                                                     std::vector<double> yy,
                                                     std::vector<double> rr){
    bool collision_detected = false;
    double dist_to_obs;
    std::vector<std::vector<double>> params(3);
    for (int k = 0; k < N_+1; k++){
        if (k != final_ind_){
            get_x(sol, k);
            for (int i = 0; i < xx.size(); i++){
                params = {{xx[i], yy[i], rr[i]}};
                dist_to_obs = std::pow(xk_[0] - xx[i], 2) + 
                              std::pow(xk_[1] - yy[i], 2);
                if (dist_to_obs < std::pow(rr[i] + 0.25, 2)){
                    nlp_->addExtraConstraint({k}, params);
                    
                    if (dist_to_obs < std::pow(rr[i]-1.0e-6, 2)){
                       collision_detected = true;
                    };
                }
            }
        }
    }
    return collision_detected;
};

bool ManyObstaclesHelper::processObstacleConstraintsCasadi(casadi::Opti& opti, 
                                        MX& xx, MX& uu, BuildingBlocks& blocks,
                                        std::vector<std::vector<double>>& sol, 
                                        std::vector<double> x, 
                                        std::vector<double> y,
                                        std::vector<double> r){
    bool collision_detected = false;
    double dist_to_obs;
    std::vector<double> params(3);
    for (int k = 0; k < N_+1; k++){
        if (k != final_ind_){
            for (int i = 0; i < x.size(); i++){
                params = {x[i], y[i], r[i]};
                dist_to_obs = std::pow(sol[0][k] - x[i], 2) + 
                              std::pow(sol[1][k] - y[i], 2);
                if (dist_to_obs < std::pow(r[i] + 0.25, 2)){
                    opti.subject_to((blocks.get_g_extra_lb(0) <= 
                        blocks.eval_g_extra({xx(Slice(), k), uu(Slice(), k), 
                                            params}, 0)[0]) <= 
                        blocks.get_g_extra_ub(0));
                    
                    if (dist_to_obs < std::pow(r[i]-1.0e-6, 2)){
                       collision_detected = true;
                    };
                }
            }
        }
    }
    return collision_detected;
};

double ManyObstaclesHelper::performAdaptiveLoop(Plotter& plotter, 
            std::vector<double>& x0, std::vector<double>& xf, int nx, int nu, 
            std::vector<double>& obstacle_x, std::vector<double>& obstacle_y, 
            std::vector<double>& obstacle_r, int max_nb_iterations, 
            std::vector<double>& timings_adaptive, 
            std::vector<std::vector<std::vector<double>>>& 
                formatted_solutions,
            std::vector<int>& nb_constraints, double& total_time,
            int& iteration_count){
    // Start with some memory allocations
    std::vector<std::vector<double>> formatted_xx;
    std::vector<std::vector<double>> formatted_uu;
    std::vector<std::optional<int>> next_ind = nlp_->getNextInd();
    std::vector<std::optional<double>> time_from_index = 
        nlp_->getTimeFromIndex();
    plotter.update(nlp_->getN(), nlp_->getFinalInd(), next_ind, 
                   time_from_index);
    std::vector<double> sol;
    iteration_count = 0;
    
    // start iterative procedure
    auto start = high_resolution_clock::now();
    for (int iteration = 0; iteration < max_nb_iterations; iteration++){
        cout<<"iteration "<<iteration<<endl;

        nb_constraints[iteration] = nlp_->getNumberOfConstraints();
        // solve the nlp
        if (iteration == 0){
            nlp_->solveNlp({{"p_g0", x0}, {"p_gT", xf}}, 
                        timings_adaptive[iteration]);
        } else {
            nlp_->solveNlp({{"p_g0", x0}, {"p_gT", xf}}, 
                        timings_adaptive[iteration], sol);
        }
        iteration_count++;
        
        // extract and format the solution
        sol = nlp_->getSolution();

        plotter.formatSolution(sol, formatted_xx, formatted_uu);
            for (int i = 0; i < nx; i++){
            formatted_solutions[iteration][i] = formatted_xx[i];
        }
        for (int i = 0; i < nu; i++){
            formatted_solutions[iteration][nx+i] = formatted_uu[i];
        }

        if (!processObstacleConstraints(sol, obstacle_x, obstacle_y, 
                                        obstacle_r)){
            break;
        };
    }
    auto stop = high_resolution_clock::now();
    total_time = double(duration_cast<nanoseconds>(stop-start).count())*1.0e-6;
    return total_time;
};

double ManyObstaclesHelper::solveCompleteNLPCasadi(BuildingBlocks& blocks, 
            Plotter& plotter, std::vector<double> x0, std::vector<double> xf, 
            int nx, int nu, int nb_steps, double dt, int nb_intervals, 
            std::vector<double>& obstacle_x, std::vector<double>& obstacle_y, 
            std::vector<double>& obstacle_r, 
            std::vector<std::vector<std::vector<double>>>& 
                formatted_solutions, double& total_time){
    int N = nb_intervals*(nb_steps-1);
    Opti opti = Opti();
    MX xx = opti.variable(nx, N+1);
    MX uu = opti.variable(nu, N);
    MX t = opti.variable(1, 1);
    std::vector<std::vector<double>> formatted_xx(nx,
        std::vector<double>(N+1));
    std::vector<std::vector<double>> formatted_uu(nu,
        std::vector<double>(N));;

    // objective
    MX obj = blocks.eval_Phi_0({xx(Slice(),0), t, 0})[0] + 
             blocks.eval_Phi_f({xx(Slice(),N), t, 0})[0];
    for (int k = 0; k < N; k++){
        obj += blocks.eval_phi({xx(Slice(), k), uu(Slice(), k), dt, t, 0})[0];
    }
    opti.minimize(obj);

    // constraints
    opti.subject_to((blocks.get_g0_lb() <= 
                    blocks.eval_g0({xx(Slice(),0), t, x0})[0]) <= 
                    blocks.get_g0_ub());
    opti.subject_to((blocks.get_gT_lb() <= 
                    blocks.eval_gT({xx(Slice(),N), xf})[0]) <= 
                    blocks.get_gT_ub());
    std::vector<double> p_obs(3);
    for (int k = 0; k < N; k++){
        opti.subject_to((blocks.get_g_fixed_lb() <= 
            blocks.eval_g_fixed({xx(Slice(), k), uu(Slice(), k), 0})[0]) <= 
            blocks.get_g_fixed_ub());
        
        for (int i = 0; i < obstacle_x.size(); i++){
            p_obs = {obstacle_x[i], obstacle_y[i], obstacle_r[i]};
            opti.subject_to((blocks.get_g_extra_lb(0) <= 
                blocks.eval_g_extra({xx(Slice(), k), uu(Slice(), k), 
                                      p_obs}, 0)[0]) <= 
                blocks.get_g_extra_ub(0));
        }
    }


    assert (nb_steps == 2 || nb_steps == 3);
    for (int j = 0; j < nb_intervals; j++){
        if (nb_steps == 2){
            MX local_xx = horzcat(xx(Slice(), j), xx(Slice(), j+1));
            opti.subject_to(
                blocks.eval_g_disc({local_xx, uu(Slice(), j), dt, t}, 
                                    nb_steps)[0] == 0);
        } else {
            MX local_xx = horzcat(xx(Slice(), 2*j), xx(Slice(), 2*j+1),
                                  xx(Slice(), 2*j+2));
            MX local_uu = horzcat(uu(Slice(), 2*j), uu(Slice(), 2*j+1));
            opti.subject_to(
                blocks.eval_g_disc({local_xx, local_uu, 
                                std::vector<double>{dt, dt}, t, 0}, 
                                nb_steps)[0] == 0);
        }
    }


    // initial guess
    for (int i = 0; i < N+1; i++){ 
        opti.set_initial(xx(Slice(), i), blocks.eval_x_init({DM(dt*i)})[0]);
    }
    for (int i = 0; i < N; i++){ 
        opti.set_initial(uu(Slice(), i), blocks.eval_u_init({DM(dt*i)})[0]);
    }
    opti.set_initial(t, blocks.get_t_init());

    // set parameters
    opti.solver("ipopt", {}, {{"print_level", 5}});

    // solve problem
    auto start = high_resolution_clock::now();
    OptiSol sol_casadi = opti.solve();
    auto stop = high_resolution_clock::now();
    DM xx_DM = sol_casadi.value(xx);
    DM uu_DM = sol_casadi.value(uu);

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
    formatted_solutions[0][0] = formatted_xx[0];
    formatted_solutions[0][1] = formatted_xx[1];
    formatted_solutions[0][2] = formatted_xx[2];
    formatted_solutions[0][3] = formatted_xx[3];
    formatted_solutions[0][4] = formatted_uu[0];
    formatted_solutions[0][5] = formatted_uu[1];
    total_time = double(duration_cast<nanoseconds>(stop-start).count())*1.0e-6;

    // plotter.plotTrajectory(formatted_solutions[0][0], formatted_solutions[0][1]);
    return total_time;
};

double ManyObstaclesHelper::performCasadiLoopReformulation(
            BuildingBlocks& blocks, Plotter& plotter, std::vector<double> x0, 
            std::vector<double> xf, int nx, int nu, int nb_steps, double dt, 
            int nb_intervals, std::vector<double>& obstacle_x, 
            std::vector<double>& obstacle_y, std::vector<double>& obstacle_r, 
            int max_nb_iterations, std::vector<double>& timings_casadi_2,
            std::vector<std::vector<std::vector<double>>>& 
                formatted_solutions,
            double& total_time, int& iteration_count){
    int N = nb_intervals*(nb_steps-1);
    iteration_count = 0;
    Opti opti = Opti();
    MX xx = opti.variable(nx, N+1);
    MX uu = opti.variable(nu, N);
    MX t = opti.variable(1, 1);
    std::vector<std::vector<double>> formatted_xx(nx,
        std::vector<double>(N+1));
    std::vector<std::vector<double>> formatted_uu(nu,
        std::vector<double>(N));;

    // objective
    MX obj = blocks.eval_Phi_0({xx(Slice(),0), t, 0})[0] + 
             blocks.eval_Phi_f({xx(Slice(),N), t, 0})[0];
    for (int k = 0; k < N; k++){
        obj += blocks.eval_phi({xx(Slice(), k), uu(Slice(), k), dt, t, 0})[0];
    }
    opti.minimize(obj);

    // constraints
    opti.subject_to((blocks.get_g0_lb() <= 
                    blocks.eval_g0({xx(Slice(),0), t, x0})[0]) <= 
                    blocks.get_g0_ub());
    opti.subject_to((blocks.get_gT_lb() <= 
                    blocks.eval_gT({xx(Slice(),N), xf})[0]) <= 
                    blocks.get_gT_ub());
    std::vector<double> p_obs(3);
    for (int k = 0; k < N; k++){
        opti.subject_to((blocks.get_g_fixed_lb() <= 
            blocks.eval_g_fixed({xx(Slice(), k), uu(Slice(), k), 0})[0]) <= 
            blocks.get_g_fixed_ub());
    }

    assert (nb_steps == 2 || nb_steps == 3);
    for (int j = 0; j < nb_intervals; j++){
        if (nb_steps == 2){
            MX local_xx = horzcat(xx(Slice(), j), xx(Slice(), j+1));
            opti.subject_to(
                blocks.eval_g_disc({local_xx, uu(Slice(), j), dt, t, 0}, 
                                    nb_steps)[0] == 0);
        } else {
            MX local_xx = horzcat(xx(Slice(), 2*j), xx(Slice(), 2*j+1),
                                  xx(Slice(), 2*j+2));
            MX local_uu = horzcat(uu(Slice(), 2*j), uu(Slice(), 2*j+1));
            opti.subject_to(
                blocks.eval_g_disc({local_xx, local_uu, 
                                std::vector<double>{dt, dt}, t, 0}, 
                                nb_steps)[0] == 0);
        }
    }

    // set parameters
    opti.solver("ipopt", {}, {{"print_level", 0}});

    // solve problem
    DM xx_DM;
    DM uu_DM;
    auto start = high_resolution_clock::now();
    for (int iteration = 0; iteration < max_nb_iterations; iteration++){
        cout<<"iteration "<<iteration<<endl;
        
        if (iteration == 0){
            for (int i = 0; i < N+1; i++){ 
                opti.set_initial(xx(Slice(), i), 
                                 blocks.eval_x_init({DM(dt*i)})[0]);
            }
            for (int i = 0; i < N; i++){ 
                opti.set_initial(uu(Slice(), i), 
                                 blocks.eval_u_init({DM(dt*i)})[0]);
            }
        } else {
            for (int k = 0; k < N-1; k++){
                opti.set_initial(xx(Slice(), k), xx_DM(Slice(), k+1));
                opti.set_initial(uu(Slice(), k), uu_DM(Slice(), k+1));
            }
            opti.set_initial(xx(Slice(), N-1), xx_DM(Slice(), N));
            opti.set_initial(xx(Slice(), N), xx_DM(Slice(), N));
            opti.set_initial(uu(Slice(), N-1), uu_DM(Slice(), N-1));
        }

        OptiSol sol_casadi = opti.solve();
        iteration_count++;
        xx_DM = sol_casadi.value(xx);
        uu_DM = sol_casadi.value(uu);
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
        formatted_solutions[iteration][0] = formatted_xx[0];
        formatted_solutions[iteration][1] = formatted_xx[1];
        formatted_solutions[iteration][2] = formatted_xx[2];
        formatted_solutions[iteration][3] = formatted_xx[3];
        formatted_solutions[iteration][4] = formatted_uu[0];
        formatted_solutions[iteration][5] = formatted_uu[1];
        timings_casadi_2[iteration] = 
            1000.0*double(sol_casadi.stats()["t_wall_total"]);

        if (!processObstacleConstraintsCasadi(opti, xx, uu, blocks, 
                                        formatted_xx, obstacle_x, obstacle_y, 
                                        obstacle_r)){
            break;
        };
    }
    auto stop = high_resolution_clock::now();
    total_time = double(duration_cast<nanoseconds>(stop-start).count())*1.0e-6;
    return total_time;
};

void ManyObstaclesHelper::get_x(std::vector<double>& sol, int k){
	int offset = (k > final_ind_)*nu_;
	for (int i = 0; i < nx_; i++){
		xk_[i] = sol[(nx_ + nu_)*k - offset + free_time_ + i];
	}
}

void ManyObstaclesHelper::get_u(std::vector<double>& sol, int k){
	int offset = (k > final_ind_)*nu_;
	for (int i = 0; i < nu_; i++){
		uk_[i] = sol[(nx_ + nu_)*k + nx_ - offset + free_time_ + i];
	}
}