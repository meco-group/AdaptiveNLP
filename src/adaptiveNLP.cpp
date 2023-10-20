#include "../include/adaptiveNLP.hpp"
#include "../include/customExceptions.hpp"
#include <cassert>
#include <chrono>

#ifndef NDEBUG
#include "../include/interfaceTester.hpp"
#endif

using namespace std;
using namespace std::chrono;

void printOptionalVector(std::vector<std::optional<double>>& v){
    for (auto e : v){
        if (e.has_value()){
            std::cout<<e.value()<<" ";
        } else {
            std::cout<<" X  ";
        }
    }
    std::cout<<endl;
}

void printOptionalVector(std::vector<std::optional<int>>& v){
    for (auto e : v){
        if (e.has_value()){
            std::cout<<std::setw(3)<<std::setfill('0')<<e.value()<<" ";
        } else {
            std::cout<<" X  ";
        }
    }
    std::cout<<endl;
}

void printSparsity(const Sparsity& sp){
    for (int i = 0; i < sp.size1(); i++){
        for (int j = 0; j < sp.size2(); j++){
            if (sp.has_nz(i,j)){
                std::cout<<"X ";
            } else {
                std::cout<<". ";
            }
        }
        std::cout<<endl;
    }
    std::cout<<endl;
}

AdaptiveNLP::AdaptiveNLP(BuildingBlocks& blocks, double T, int Nmax, 
                         int max_nb_extra_instances_)
                        :Nmax_(Nmax), nx_(blocks.nx_), nu_(blocks.nu_),
                         T_(T), free_time_(blocks.free_time_),
                         blocks_(&blocks),
                         bookkeeper_(blocks, Nmax_, 
                                    getMaxNbConstraints(
                                            max_nb_extra_instances_), 
                                    max_nb_extra_instances_){
    assert (blocks.check_complete());
    N_ = 0;

    ipopt_app_ = new IpoptApplication();
    interface_ = new NLPInterface(N_, bookkeeper_, blocks);
    #ifndef NDEBUG
    interfaceTester_ = InterfaceTester(*interface_);
    #endif

    k_vals_ = std::vector<int>(blocks_->max_nb_steps_);
    
	nz_rows_ = std::vector<int>(blocks_->max_nnz_);
	nz_cols_ = std::vector<int>(blocks_->max_nnz_);
    nz_indices_ = std::vector<int>(blocks_->max_nnz_);

    selected_rows_nxnu_ = std::vector<casadi_int>(nx_ + nu_);
    selected_rows_nx_ = std::vector<casadi_int>(nx_);
    selected_cols_nxnu_ = std::vector<casadi_int>(nx_ + nu_);
    selected_cols_nx_ = std::vector<casadi_int>(nx_);
};

void AdaptiveNLP::initTimeSteps(std::vector<int> nks, std::vector<double> tt){
    // check input
    if (N_ != 0){ throw illegalChangeToNLP("Cannot call initTimeSteps if "
        "some time-steps have already been initialized.");}
    if (nks.size() == 0){ throw illegalChangeToNLP("Invalid input to "
        "function to change NLP (nks.size() == 0)");}
    for (int i = 0; i < nks.size(); i++){
        if (!blocks_->nb_steps_mapping_[nks[i]].has_value()){
            throw illegalChangeToNLP("Invalid input to function to change "
                "NLP (no discretization function exists with required "
                "number of time-steps)");
        }
        if (nks[i] <= 1){
            throw illegalChangeToNLP("Invalid input to function to change "
                "NLP (nks[i] <= 1)");
        };
    }

    // compute total number of time-steps based on input
    int  nb_to_add = 1;
    for (int i = 0; i < nks.size(); i++){
        nb_to_add += nks[i]-1;
    }
    
    // check provided timings
    if (tt.size() > 0 && tt.size() != nb_to_add){
        throw illegalChangeToNLP("Invalid time-vector provided. The size of "
                                 "the time vector should equal 1+sum(i-1) for "
                                 "i in nks");
    }
    if (tt.size() > 0){
        for (int i = 0; i < tt.size()-1; i++){
            if (tt[i] >= tt[i+1]){
                throw illegalChangeToNLP("Invalid time-vector provided. The "
                                        "times should be monotonically "
                                        "increasing.");
            }
        }
    }

    if (Nmax_ < nb_to_add-1){
        throw illegalChangeToNLP(
            "new number of primal variables exceeds allocated space");
    }

    // update bookkeeping information
    for (int k = 0; k < nb_to_add-1; k++){
        if (tt.size() == 0){
            bookkeeper_.time_from_ind_[k] = T_/(nb_to_add-1)*k;
        } else { bookkeeper_.time_from_ind_[k] = tt[k];}
        bookkeeper_.next_ind_[k] = k+1;
    }
    if (tt.size() == 0){
        bookkeeper_.time_from_ind_[nb_to_add-1] = T_;
    } else { bookkeeper_.time_from_ind_[nb_to_add-1] = tt[nb_to_add-1];}
    bookkeeper_.final_ind_ = nb_to_add - 1;
    bookkeeper_.final_ind_old_ = bookkeeper_.final_ind_;

    ///////////////////////////////////
    // hessian of objective function //
    ///////////////////////////////////
    addCompleteHessianContribution(blocks_->Phi_0_H_.sparsity_out(0), {0}, 
                                   free_time_, bookkeeper_.hess_nz_ind_Phi_0_);
    addCompleteHessianContribution(blocks_->Phi_f_H_.sparsity_out(0),
                                   {bookkeeper_.final_ind_}, free_time_, 
                                   bookkeeper_.hess_nz_ind_Phi_f_);

    ////////////////////////
    // initial constraint //
    ////////////////////////
    // lower and upper bounds on constraint value
    std::copy(blocks_->g0_lb_.begin(), blocks_->g0_lb_.end(),
              bookkeeper_.lb_container_.begin());
    std::copy(blocks_->g0_ub_.begin(), blocks_->g0_ub_.end(),
              bookkeeper_.ub_container_.begin());
    bookkeeper_.next_free_ind_ += blocks_->nb_g0_;
    // jacobian sparsity
    addCompleteJacobianContribution(blocks_->g0_J_.sparsity_out(0), {0},
                                    free_time_, 0, bookkeeper_.jac_nz_ind_g0_);
    // hessian sparsity
    addCompleteHessianContribution(blocks_->g0_H_.sparsity_out(0), {0},
                                   free_time_, bookkeeper_.hess_nz_ind_g0_);

    //////////////////////
    // final constraint //
    //////////////////////
    // lower and upper bound on constraint value
    std::copy(blocks_->gT_lb_.begin(), blocks_->gT_lb_.end(),
              bookkeeper_.lb_container_.begin() + bookkeeper_.next_free_ind_);
    std::copy(blocks_->gT_ub_.begin(), blocks_->gT_ub_.end(),
              bookkeeper_.ub_container_.begin() + bookkeeper_.next_free_ind_);
    bookkeeper_.next_free_ind_ += blocks_->nb_gT_;
    // jacobian sparsity
    addCompleteJacobianContribution(blocks_->gT_J_.sparsity_out(0),
                                    {bookkeeper_.final_ind_}, false,
                                    blocks_->nb_g0_, 
                                    bookkeeper_.jac_nz_ind_gT_);
    // hessian sparsity
    addCompleteHessianContribution(blocks_->gT_H_.sparsity_out(0),
                                   {bookkeeper_.final_ind_}, false,
                                   bookkeeper_.hess_nz_ind_gT_);

    // stage contributions
    addStageContributions(0, nb_to_add-1);

    // discretization constraints
    addDiscretizationContributions(0, nks);

    N_ += nb_to_add - 1;
    interface_->updateN(N_);
};

void AdaptiveNLP::addTimeSteps(std::vector<int> nks, std::vector<double> tt){
    // check input
    int nb_to_add = checkInputToNLPUpdate(nks);

    // check provided timings
    if (tt.size() > 0 && tt.size() != nb_to_add){
        throw illegalChangeToNLP("Invalid time-vector provided. The size of "
                                 "the time vector should equal sum(i-1) for "
                                 "i in nks");
    }
    if (tt.size() > 0){
        for (int i = 0; i < tt.size()-1; i++){
            if (tt[i] >= tt[i+1]){
                throw illegalChangeToNLP("Invalid time-vector provided. The "
                                        "times should be monotonically "
                                        "increasing.");
            }
        }
    }

    // remove the final ind from the current chain of time-steps
    // get the second-to-last index
    int almost_final_index = 0;
    while (bookkeeper_.next(almost_final_index) != bookkeeper_.final_ind_){
        almost_final_index = bookkeeper_.next(almost_final_index);
    }
   
    // change the next-value of the second-to-last index
    bookkeeper_.next_ind_[almost_final_index] = N_+1;

    // update bookkeeping information
    double new_dt = T_/(N_+nb_to_add);

    for (int k = 1; k < nb_to_add; k++){
        bookkeeper_.next_ind_[N_+k] = N_+k+1;
    }
    bookkeeper_.next_ind_[N_+nb_to_add] = bookkeeper_.final_ind_;

    // update all time-steps
    int k = 0;
    int chronological_k = 0;
    bookkeeper_.time_from_ind_[k] = chronological_k*new_dt;
    while (bookkeeper_.next_ind_[k].has_value()){
        k = bookkeeper_.next_ind_[k].value();
        chronological_k++;
        bookkeeper_.time_from_ind_[k] = chronological_k*new_dt;
    }

    // overwrite the discretization constraint
    int k_init = bookkeeper_.ind_x_d_[almost_final_index].value();
    getIndices(k_init, bookkeeper_.nb_g_d_[k_init].value());
    addCompleteJacobianContribution(blocks_->g_disc_J_[
                                        bookkeeper_.ind_disc_block_[k_init]
                                        .value()].sparsity_out(0),
                                    std::vector<int>(k_vals_.begin(), 
                                        k_vals_.begin() + 
                                        bookkeeper_.nb_g_d_[k_init].value()), 
                                    free_time_,
                                    bookkeeper_.ind_g_d_[k_init].value(), 
                                    bookkeeper_.jac_nz_ind_g_disc_[k_init]);

    addCompleteHessianContribution(blocks_->g_disc_H_[
                                        bookkeeper_.ind_disc_block_[k_init]
                                        .value()].sparsity_out(0),
                                   std::vector<int>(k_vals_.begin(), 
                                        k_vals_.begin() + 
                                        bookkeeper_.nb_g_d_[k_init].value()), 
                                   free_time_,
                                   bookkeeper_.hess_nz_ind_g_disc_[k_init]);

    // add stage-cost and fixed constraints
    addStageContributions(N_+1, nb_to_add);
   
    // dynamical constraints
    addDiscretizationContributions(N_+1, nks);

    N_ += nb_to_add;
    interface_->updateN(N_);

    clearStructuralZeros();
};

void AdaptiveNLP::appendTimeSteps(std::vector<int> nks, 
                                  std::vector<double> tt){
    // Note on implementation:
    //  - The final new sequence of time-steps is inserted in between
    //    the second-to-last timestep and the final_index-time-step. Hence,
    //    there is no need to update the final constraint or final cost.
    //  - Then, the time of the final index is updated accordingly.

    // check input
    int nb_to_add = checkInputToNLPUpdate(nks);
    
    // check provided timings
    if (tt.size() > 0 && tt.size() != nb_to_add){
        throw illegalChangeToNLP("Invalid time-vector provided. The size of "
                                 "the time vector should equal sum(i-1) for "
                                 "i in nks");
    }
    if (tt.size() > 0){
        for (int i = 0; i < tt.size()-1; i++){
            if (tt[i] >= tt[i+1]){
                throw illegalChangeToNLP("Invalid time-vector provided. The "
                                        "times should be monotonically "
                                        "increasing.");
            }
        }
    }

    // remove the final ind from the current chain of time-steps
    // get the second-to-last index
    int almost_final_index = 0;
    while (bookkeeper_.next(almost_final_index) != bookkeeper_.final_ind_){
        almost_final_index = bookkeeper_.next(almost_final_index);
    }
   
    // change the next-value of the second-to-last index (such that it points 
    // to the appending sequence)
    bookkeeper_.next_ind_[almost_final_index] = N_+1;

    // update bookkeeping information
    double dt = T_/N_;
    T_ += nb_to_add*dt;

    for (int k = 1; k < nb_to_add; k++){
        if (tt.size() == 0){
            bookkeeper_.time_from_ind_[N_+k] = 
                bookkeeper_.time_from_ind_[almost_final_index].value() + dt*k;
        } else { bookkeeper_.time_from_ind_[N_+k] = tt[k-1];}
        bookkeeper_.next_ind_[N_+k] = N_+k+1;
    }
    bookkeeper_.next_ind_[N_+nb_to_add] = bookkeeper_.final_ind_;
    bookkeeper_.time_from_ind_[N_+nb_to_add] = 
        bookkeeper_.time_from_ind_[almost_final_index].value() + 
        dt*nb_to_add;
    bookkeeper_.time_from_ind_[bookkeeper_.final_ind_] = 
        bookkeeper_.time_from_ind_[bookkeeper_.final_ind_].value() + 
        dt*nb_to_add;

    // overwrite the discretization constraint
    int k_init = bookkeeper_.ind_x_d_[almost_final_index].value();
    getIndices(k_init, bookkeeper_.nb_g_d_[k_init].value());
    addCompleteJacobianContribution(blocks_->g_disc_J_[
                                        bookkeeper_.ind_disc_block_[k_init]
                                        .value()].sparsity_out(0),
                                    std::vector<int>(k_vals_.begin(), 
                                        k_vals_.begin() + 
                                        bookkeeper_.nb_g_d_[k_init].value()), 
                                    free_time_,
                                    bookkeeper_.ind_g_d_[k_init].value(), 
                                    bookkeeper_.jac_nz_ind_g_disc_[k_init]);

    addCompleteHessianContribution(blocks_->g_disc_H_[
                                        bookkeeper_.ind_disc_block_[k_init]
                                        .value()].sparsity_out(0),
                                   std::vector<int>(k_vals_.begin(), 
                                        k_vals_.begin() + 
                                        bookkeeper_.nb_g_d_[k_init].value()), 
                                   free_time_,
                                   bookkeeper_.hess_nz_ind_g_disc_[k_init]);

    // add stage-cost and fixed constraints
    addStageContributions(N_+1, nb_to_add);
   
    // dynamical constraints
    addDiscretizationContributions(N_+1, nks);

    N_ += nb_to_add;

    interface_->updateN(N_);
    clearStructuralZeros();
};

void AdaptiveNLP::reduceHorizon(int nb_intervals){

    // compute new final index and count number of time steps to be removed
    int nb_steps_counter = 0;
    int new_final_ind = bookkeeper_.final_ind_;
    for (int i = 0; i < nb_intervals; i++){
        new_final_ind = bookkeeper_.ind_x_d_[
            bookkeeper_.previous(new_final_ind)].value();
        nb_steps_counter += bookkeeper_.nb_g_d_[new_final_ind].value() - 1;
    }
    bookkeeper_.final_ind_old_ = bookkeeper_.final_ind_;
    bookkeeper_.final_ind_ = new_final_ind;
    bookkeeper_.next_ind_[new_final_ind].reset();
    N_ -= nb_steps_counter;
    interface_->updateN(N_);

    // update final objective function contribution
    addCompleteHessianContribution(blocks_->Phi_f_H_.sparsity_out(0),
                                   {bookkeeper_.final_ind_}, free_time_, 
                                   bookkeeper_.hess_nz_ind_Phi_f_);

    // update final state constraint
    // jacobian sparsity
    addCompleteJacobianContribution(blocks_->gT_J_.sparsity_out(0),
                                    {bookkeeper_.final_ind_}, false,
                                    blocks_->nb_g0_, 
                                    bookkeeper_.jac_nz_ind_gT_);
    // hessian sparsity
    addCompleteHessianContribution(blocks_->gT_H_.sparsity_out(0),
                                   {bookkeeper_.final_ind_}, false,
                                   bookkeeper_.hess_nz_ind_gT_);
    
    bookkeeper_.remove_traces_old_vars_ = true;
};

void AdaptiveNLP::changeIntervalDiscretization(int k, std::vector<int> nks,
                                               std::vector<double> tt){
    int k_init = k;
    // check input
    if (k < 0 || k > N_ || k != bookkeeper_.ind_x_d_[k].value()){
        throw illegalChangeToNLP("cannot change discretization constraints "
                                 "of interval at invalid k-index (k = " + 
                                 std::to_string(k) + ")");
    }
    if (nks.size() == 1 && nks[0] == bookkeeper_.nb_g_d_[k].value()){
        return;
    }
    int nb_to_add_bruto = checkInputToNLPUpdate(nks) - 1;
    
    // check provided timings
    if (tt.size() > 0 && tt.size() != nb_to_add_bruto){
        throw illegalChangeToNLP("Invalid time-vector provided. The size of "
                                 "the time vector should equal (sum(i-1) for "
                                 "i in nks) - 1");
    }
    if (tt.size() > 0){
        for (int i = 0; i < tt.size()-1; i++){
            if (tt[i] >= tt[i+1]){
                cout<<tt<<endl;
                throw illegalChangeToNLP("Invalid time-vector provided. The "
                                        "times should be monotonically "
                                        "increasing.");
            }
        }
    }

    // compute the total change to the number of time-steps
    int nb_to_add_netto = nb_to_add_bruto - bookkeeper_.nb_g_d_[k].value() + 2;

    // get the index of the first time-step of the next interval (in the old
    // discretization)
    int k_attach = k;
    for (int i = 0; i < bookkeeper_.nb_g_d_[k].value()-1; i++){
        k_attach = bookkeeper_.next(k_attach);
    }

    if (tt[tt.size()-1] >= bookkeeper_.time_from_ind_[k_attach].value()){
        throw illegalChangeToNLP("Invalid time-vector provided. The last "
                                 "element of the time vector should be lower "
                                 "than the time of the index to which the "
                                 "new discretization is attached");
    }
    if (tt[0] <= bookkeeper_.time_from_ind_[k_init].value()){
        throw illegalChangeToNLP("Invalid time-vector provided. The first "
                                 "element of the time vector should be higher "
                                 "than the time of the index at which the "
                                 "new discretization is started");
    }

    // clear some bookkeeping
    for (auto& e : bookkeeper_.jac_nz_ind_g_disc_[k_init]){e.reset();}
    for (auto& e : bookkeeper_.hess_nz_ind_g_disc_[k_init]){e.reset();}

    // case where there will be less time-steps
    if (nb_to_add_netto < 0){
        bookkeeper_.remove_traces_old_vars_ = true;

        // make the last kept time-step point to the first time-step of the
        // next interval in the old discretization
        k = k_init;
        for (int i = 0; i < nb_to_add_bruto; i++){
            k = bookkeeper_.next(k);
            if (tt.size() > 0){bookkeeper_.time_from_ind_[k] = tt[i];}
        }
        bookkeeper_.next_ind_[k] = k_attach;

        // update time_from_ind
        if (tt.size() == 0){
            double dt = (bookkeeper_.time_from_ind_[k_attach].value() -
                    bookkeeper_.time_from_ind_[k_init].value())/
                    (nb_to_add_bruto+1);
            k = k_init;
            double time = bookkeeper_.time_from_ind_[k_init].value();
            int counter = 0;
            while (k != k_attach){
                bookkeeper_.time_from_ind_[k] = time + counter*dt;
                k = bookkeeper_.next(k);
                counter++;
            }
        } else {
            k = bookkeeper_.next(k_init);
            int counter = 0;
            while (k!= k_attach){
                bookkeeper_.time_from_ind_[k] = tt[counter];
                k = bookkeeper_.next(k);
                counter++;
            }
        }

        // dynamical constraints
        addDiscretizationContributions(k_init, nks);

    // case where time-steps have to be added
    } else if (nb_to_add_netto > 0){
        // add new time-steps

        int almost_final_index = 0;
        while (bookkeeper_.next(almost_final_index) != k_attach){
            almost_final_index = bookkeeper_.next(almost_final_index);
        }
    
        // change the next-value of the second-to-last index (such that it points 
        // to the appending sequence)
        bookkeeper_.next_ind_[almost_final_index] = N_+1;

        // update next-chain
        for (int k = 1; k < nb_to_add_netto; k++){
            bookkeeper_.next_ind_[N_+k] = N_+k+1;
        }
        bookkeeper_.next_ind_[N_+nb_to_add_netto] = k_attach;

        // update time_from_ind
        double dt = (bookkeeper_.time_from_ind_[k_attach].value() -
                    bookkeeper_.time_from_ind_[k_init].value())/
                        (nb_to_add_bruto+1);
        int k = bookkeeper_.next(k_init);
        double time = bookkeeper_.time_from_ind_[k_init].value();
        int counter = 1;
        while (k != k_attach){
            if (tt.size() == 0){
                bookkeeper_.time_from_ind_[k] = time + counter*dt;
            } else {bookkeeper_.time_from_ind_[k] = tt[counter-1];}
            k = bookkeeper_.next(k);
            counter++;
        }

        // add stage-cost and fixed constraints
        addStageContributions(N_+1, nb_to_add_netto);
    
        // dynamical constraints
        addDiscretizationContributions(k_init, nks);
    
    // case where the number of time-steps stays the same
    } else {
        // dynamical constraints
        addDiscretizationContributions(k_init, nks);
    }

    N_ += nb_to_add_netto;
    interface_->updateN(N_);
};

void AdaptiveNLP::addExtraConstraint(std::vector<int> inds, 
                                     std::vector<int> constraint_ind,
                                     std::vector<int> instance_ind,
                                     std::vector<std::vector<double>> 
                                                                parameters){
    // check input
    if (constraint_ind.size() > 1 && constraint_ind.size() != inds.size() ||
            parameters.size() > 1 && parameters.size() != inds.size() ||
            instance_ind.size() > 1 && instance_ind.size() != inds.size() ||
            constraint_ind.size() == 1 && parameters.size() != 1 ||
            constraint_ind.size() != 1 && parameters.size() == 1 ||
            constraint_ind.size() == 1 && instance_ind.size() != 1 ||
            constraint_ind.size() != 1 && instance_ind.size() == 1){
        throw illegalChangeToNLP("invalid lengths provided for constraint "
                                 "indices and parameters in "
                                 "'addExtraConstraint()'");
    }

    bool use_same_constraint_ind = (constraint_ind.size() == 1);
    bool reset_instance_ind = false;
    if (use_same_constraint_ind && instance_ind[0] == -1){
        reset_instance_ind = true;
    }
    int ind_ptr = 0;

    for (int i = 0; i < inds.size(); i++){
        if (!use_same_constraint_ind){ ind_ptr = i;};

        if (inds[i] < 0 || inds[i] > N_+1 || 
                inds[i] == bookkeeper_.final_ind_){
            throw illegalChangeToNLP("Cannot add extra constraint to invalid "
                                     "k-index (k = " + 
                                    std::to_string(inds[i]) + ")");
        }
        if (!(blocks_->nb_p_g_extra_[constraint_ind[ind_ptr]] == 0 ||
                    parameters[ind_ptr].size() == 
                        blocks_->nb_p_g_extra_[constraint_ind[ind_ptr]])){
            throw illegalChangeToNLP("Invalid parameter sizes for extra "
                                     "constraints to be added.");
        };
    }

    // add constraints
    ind_ptr = 0;
    for (int k = 0; k < inds.size(); k++){
        if (!use_same_constraint_ind){ ind_ptr = k;};
        if (reset_instance_ind){ instance_ind[ind_ptr] = -1;}

        // check if there is a free slot for another instance of this 
        // constraint
        if (instance_ind[ind_ptr] == -1 && 
                !getFreeInstanceSlot(inds[k], constraint_ind[ind_ptr], 
                                     instance_ind[ind_ptr])){
            throw illegalChangeToNLP("Cannot add another instance of "
                                     "constraint with index " + 
                                     std::to_string(
                                     constraint_ind[ind_ptr]) +
                                     " to time-index k = " +
                                     std::to_string(inds[k]));
        }

        if (instance_ind[ind_ptr] != -1 && 
                bookkeeper_.ind_g_e_.count({inds[k], constraint_ind[ind_ptr], 
                                      instance_ind[ind_ptr]}) > 0){
            throw illegalChangeToNLP("Cannot add instance (" + 
                                     std::to_string(instance_ind[ind_ptr]) + 
                                     ") of extra constraint (" +
                                     std::to_string(constraint_ind[ind_ptr]) +
                                     ") to unavailable slot at time-step " +
                                     std::to_string(inds[k]));
        }

        // add parameter to interface
        if (blocks_->nb_p_g_extra_[constraint_ind[ind_ptr]] > 0){
            std::copy(parameters[ind_ptr].begin(), parameters[ind_ptr].end(),
                interface_->
                    p_g_extra_[inds[k]][constraint_ind[ind_ptr]]
                        [instance_ind[ind_ptr]].begin());
        }
        
        // update constraint bounds
        std::copy(blocks_->g_extra_lb_[constraint_ind[ind_ptr]].begin(),
                  blocks_->g_extra_lb_[constraint_ind[ind_ptr]].end(),
                  bookkeeper_.lb_container_.begin() + 
                    bookkeeper_.next_free_ind_);
        std::copy(blocks_->g_extra_ub_[constraint_ind[ind_ptr]].begin(),
                  blocks_->g_extra_ub_[constraint_ind[ind_ptr]].end(),
                  bookkeeper_.ub_container_.begin() + 
                    bookkeeper_.next_free_ind_);
        
        // constraint bookkeeping
        bookkeeper_.ind_g_e_[{inds[k], constraint_ind[ind_ptr], 
                              instance_ind[ind_ptr]}] = 
                                    bookkeeper_.next_free_ind_;
        bookkeeper_.next_free_ind_ += blocks_->
                                        nb_g_extra_[constraint_ind[ind_ptr]];
        bookkeeper_.nb_g_extra_applied_[constraint_ind[ind_ptr]] += 
            blocks_->nb_g_extra_[constraint_ind[ind_ptr]];

        // jacobian
        addCompleteJacobianContribution(
            blocks_->g_extra_J_[constraint_ind[ind_ptr]].sparsity_out(0), 
            {inds[k]}, false,
            bookkeeper_.ind_g_e_[{inds[k], constraint_ind[ind_ptr], 
                                  instance_ind[ind_ptr]}],
            bookkeeper_.
                jac_nz_ind_g_extra_[inds[k]][constraint_ind[ind_ptr]]
                [instance_ind[ind_ptr]]);

        // hessian
        addCompleteHessianContribution(
            blocks_->g_extra_H_[constraint_ind[ind_ptr]].sparsity_out(0), 
            {inds[k]}, false, bookkeeper_.
                hess_nz_ind_g_extra_[inds[k]][constraint_ind[ind_ptr]]
                    [instance_ind[ind_ptr]]);
       }
};

void AdaptiveNLP::removeExtraConstraint(std::vector<int> inds, 
                                        std::vector<int> constraint_ind, 
                                        std::vector<int> instance_ind){
    assert (inds.size() == constraint_ind.size());
    assert (inds.size() == instance_ind.size());

    for (int k = 0; k < inds.size(); k++){
        // if there is no valid constraint index or instance index, try to
        // retrieve them. If this fails, throw an exception.
        // Also throw an expection if there is no constraint for the given
        // indices
        if (( constraint_ind[k] == -1 && 
               !getFirstConstraintAndInstanceInd(inds[k], constraint_ind[k], 
                                                 instance_ind[k]) )
            || ( instance_ind[k] == -1 && 
                !getFirstInstanceInd(inds[k], constraint_ind[k], 
                                     instance_ind[k]) )
            || ( bookkeeper_.ind_g_e_.count({inds[k], constraint_ind[k], 
                                             instance_ind[k]}) == 0)){
            // throw illegalConstraintRemovalException();
            std::cout<<"skipping removal of a constraint"<<endl;
            continue;
        }

        // if there is a constraint, remove it
        bookkeeper_.ind_g_e_.erase({inds[k], constraint_ind[k], 
                                    instance_ind[k]});
        bookkeeper_.nb_g_extra_applied_[constraint_ind[k]] -=
            blocks_->nb_g_extra_[constraint_ind[k]];
    }
}

void AdaptiveNLP::removeExtraConstraint(std::vector<int> inds, 
                                        int constraint_ind, 
                                        int instance_ind){   
    for (int k = 0; k < inds.size(); k++){
        // if there is no valid constraint index or instance index, try to
        // retrieve them. If this fails, throw an exception.
        // Also throw an expection if there is no constraint for the given
        // indices
        if (inds[k] == bookkeeper_.final_ind_){
            throw illegalConstraintRemovalException();
        }

        if (( constraint_ind == -1 && 
               !getFirstConstraintAndInstanceInd(inds[k], constraint_ind, 
                                                 instance_ind) )
            || ( instance_ind == -1 && 
                !getFirstInstanceInd(inds[k], constraint_ind, 
                                     instance_ind) )
            || ( bookkeeper_.ind_g_e_.count({inds[k], constraint_ind, 
                                             instance_ind}) == 0)){
            continue;
        }

        // if there is a constraint, remove it
        bookkeeper_.ind_g_e_.erase({inds[k], constraint_ind, instance_ind});
        bookkeeper_.nb_g_extra_applied_[constraint_ind] -= 
            blocks_->nb_g_extra_[constraint_ind];
    }
}

int AdaptiveNLP::solveNlp(std::map<std::string,
                          std::vector<double>> parameters,
                          double& computation_time){
    assert (N_ > 0);
    if (bookkeeper_.remove_traces_old_constraints_ || 
            bookkeeper_.remove_traces_old_vars_){
        throw new illegalCallToSolve();
    }
    
    interface_->updateParameters(parameters);
    if (!interface_->checkParameters()){
        std::cerr<<"ERROR: Invalid parameters provided!"<<endl;
        std::cerr<<"Make sure to provide the correct number of parameters: "
            <<endl;
        std::cerr << "Phi_0:   " << blocks_->nb_p_Phi_0_ << std::endl;
        std::cerr << "Phi_f:   " << blocks_->nb_p_Phi_f_ << std::endl;
        std::cerr << "phi:     " << blocks_->nb_p_phi_ << std::endl;
        std::cerr << "g0:      " << blocks_->nb_p_g0_ << std::endl;
        std::cerr << "gT:      " << blocks_->nb_p_gT_ << std::endl;
        std::cerr << "g_fixed: " << blocks_->nb_p_g_fixed_ << std::endl;
        std::cerr << "g_disc:  " << blocks_->nb_p_g_disc_ << std::endl;
        
        throw illegalParameterProvided();
    }
    ipopt_app_->Options()->SetIntegerValue("print_level", print_level_);
    // ipopt_app_->Options()->SetNumericValue("tol", 1.0e-15);

    auto start = high_resolution_clock::now();

    ipopt_app_->Initialize();
    ApplicationReturnStatus status = ipopt_app_->OptimizeTNLP(interface_);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    computation_time = double(duration.count())/1000.0;

    return (int) status;
};

std::vector<double> AdaptiveNLP::getSolution(){
    std::vector<double> res(interface_->sol_x_.size());
    std::copy(interface_->sol_x_.begin(), interface_->sol_x_.end(),
              res.begin());
    return res;
}

#ifndef NDEBUG
void AdaptiveNLP::testInterface(std::vector<double>& x){    
    // std::cout<<"[AdaptiveNLP::testInterface()] evaluating f"<<endl;
    // interfaceTester_.eval_f(x);
    
    // std::cout<<"[AdaptiveNLP::testInterface()] evaluating grad_f"<<endl;
    // interfaceTester_.eval_grad_f(x);

    std::cout<<"[AdaptiveNLP::testInterface()] evaluating constraints"<<endl;
    interfaceTester_.eval_g(x);

    std::cout<<"[AdaptiveNLP::testInterface()] evaluating jacobian"<<endl;
    interfaceTester_.eval_jac_g(x);

    // std::cout<<"[AdaptiveNLP::testInterface()] evaluating hessian"<<endl;
    // std::vector<double> lambda = std::vector<double>(bookkeeper_.next_free_ind_,
    //                                                  1.0);
    // interfaceTester_.eval_h(x, lambda);
};
#endif

int AdaptiveNLP::checkInputToNLPUpdate(std::vector<int>& nks){
    if (N_ == 0){
        throw illegalChangeToNLP("Cannot add time-steps without having "
                                 "initialized some time-steps. Please call "
                                 "initTimeSteps first.");
    }
    if (nks.size() == 0){ throw illegalChangeToNLP("Invalid input to "
                                                   "function to change NLP "
                                                   "(nks.size() == 0)");}
    for (int i = 0; i < nks.size(); i++){
        if (nks[i] <= 1){ throw illegalChangeToNLP("Invalid input to "
                                                   "function to change NLP "
                                                   "(nks[i] <= 1)");}
        if (!blocks_->nb_steps_mapping_[nks[i]].has_value()){
            throw illegalChangeToNLP("Invalid input to function to change "
                "NLP (no discretization function exists with required "
                "number of time-steps)");
        }
    }

    // compute total number of time-steps to add based on input
    int  nb_to_add = 0;
    for (int i = 0; i < nks.size(); i++){
        nb_to_add += nks[i]-1;
    }
    if (Nmax_ < N_ + nb_to_add){
        throw illegalChangeToNLP(
            "new number of primal variables exceeds allocated space");
    }
    return nb_to_add;
}

void AdaptiveNLP::addStageContributions(int k_start, int nb_of_stages){
    for (int k = k_start; k < k_start+nb_of_stages; k++){
        // hessian constribution to stage cost 
        addCompleteHessianContribution(blocks_->phi_H_.sparsity_out(0), {k},
                                       free_time_,
                                       bookkeeper_.hess_nz_ind_phi_[k]);

        // lower and upper bound on fixed constraint value
        std::copy(blocks_->g_fixed_lb_.begin(), blocks_->g_fixed_lb_.end(),
                  bookkeeper_.lb_container_.begin() + 
                    bookkeeper_.next_free_ind_);
        std::copy(blocks_->g_fixed_ub_.begin(), blocks_->g_fixed_ub_.end(),
                  bookkeeper_.ub_container_.begin() + 
                    bookkeeper_.next_free_ind_);

        // bookkeeping
        bookkeeper_.ind_g_f_[k] = bookkeeper_.next_free_ind_;
        bookkeeper_.next_free_ind_ += blocks_->nb_g_fixed_;
       
        // jacobian sparsity
        addCompleteJacobianContribution(blocks_->g_fixed_J_.sparsity_out(0),
                                        {k}, false, 
                                        bookkeeper_.ind_g_f_[k].value(),
                                        bookkeeper_.jac_nz_ind_g_fixed_[k]);

        // hessian sparsity
        addCompleteHessianContribution(blocks_->g_fixed_H_.sparsity_out(0),
                                        {k}, false,
                                        bookkeeper_.hess_nz_ind_g_fixed_[k]);
    }
};

void AdaptiveNLP::addDiscretizationContributions(int k_start, 
                                                 std::vector<int>& nks){
    int k = k_start;
    for (int i = 0; i < nks.size(); i++){
        // lower and upper bound on constraint value
        std::fill(bookkeeper_.lb_container_.begin() + 
                    bookkeeper_.next_free_ind_,
                  bookkeeper_.lb_container_.begin() + 
                    bookkeeper_.next_free_ind_ + 
                    blocks_->nb_g_disc_[blocks_->
                        nb_steps_mapping_[nks[i]].value()], 0);
        std::fill(bookkeeper_.ub_container_.begin() + 
                    bookkeeper_.next_free_ind_,
                  bookkeeper_.ub_container_.begin() + 
                    bookkeeper_.next_free_ind_ + 
                    blocks_->nb_g_disc_[blocks_->
                        nb_steps_mapping_[nks[i]].value()], 0);

        // bookkeeping
        int k_local = k;
        for (int j = 0; j < nks[i]-1; j++){
            bookkeeper_.ind_g_d_[k_local] = bookkeeper_.next_free_ind_;
            bookkeeper_.nb_g_d_[k_local] = nks[i];
            bookkeeper_.ind_disc_block_[k_local] = 
                blocks_->nb_steps_mapping_[nks[i]];
            bookkeeper_.ind_x_d_[k_local] = k;
            k_local = bookkeeper_.next(k_local);
        }
        getIndices(k, nks[i]);
        bookkeeper_.next_free_ind_ += 
            blocks_->nb_g_disc_[blocks_->nb_steps_mapping_[nks[i]].value()];

        // jacobian sparsity
        addCompleteJacobianContribution(blocks_->g_disc_J_[blocks_->
                                            nb_steps_mapping_[nks[i]].value()]
                                            .sparsity_out(0),
                                        std::vector<int>(k_vals_.begin(), 
                                                         k_vals_.begin() + 
                                                         nks[i]), 
                                        free_time_,
                                        bookkeeper_.ind_g_d_[k].value(), 
                                        bookkeeper_.jac_nz_ind_g_disc_[k]);

        // hessian sparsity
        addCompleteHessianContribution(blocks_->g_disc_H_[blocks_->
                                            nb_steps_mapping_[nks[i]].value()]
                                            .sparsity_out(0),
                                       std::vector<int>(k_vals_.begin(), 
                                                        k_vals_.begin() + 
                                                        nks[i]), 
                                       free_time_,
                                       bookkeeper_.hess_nz_ind_g_disc_[k]);
        k = k_vals_[nks[i]-1];
    }
};

void AdaptiveNLP::getIndices(int k_init, int nk){
    k_vals_[0] = k_init;
    for (int i = 1; i < nk; i++){
        k_init = bookkeeper_.next(k_init);
        k_vals_[i] = k_init;
    }
};

int AdaptiveNLP::getNonzeroIndices(const Sparsity& sp,
                                   bool only_lower_triangular, 
                                   int row_offset, int col_offset){
	std::vector<casadi_int> colind = sp.get_colind();
	std::vector<casadi_int> rowind = sp.get_row();
	int idx_ptr = 0;
	for (int col_ptr = 1; col_ptr < colind.size(); col_ptr++){
		for (int row_ptr = colind[col_ptr-1]; row_ptr < colind[col_ptr]; 
				row_ptr++){
            if (!only_lower_triangular || 
                    row_offset + rowind[row_ptr] >= col_offset + col_ptr-1){
                nz_cols_[idx_ptr] = col_ptr-1;
                nz_rows_[idx_ptr] = rowind[row_ptr];
                idx_ptr++;
            }
		}
	}
    return idx_ptr;
};

int AdaptiveNLP::addElementaryContribution(Sparsity& sp, int row_offset,
                                    int col_offset, bool hessian,
                                    std::vector<std::optional<int>>& row_ind,
                                    std::vector<std::optional<int>>& col_ind){
    int nnz;
    if (hessian){
        nnz = bookkeeper_.hess_nnz_;
    } else {
        nnz = bookkeeper_.jac_nnz_;
    };

    // get the row and column indices of the nonzero elements
    int nnz_added = getNonzeroIndices(sp, hessian, row_offset, col_offset);

    // loop over all nonzers to be added
    bool new_nz;
    int ind;
    for (int j = 0; j < nnz_added; j++){
        new_nz = true;
        ind = 0;

        // if there exists already a nonzero at the same location, find its
        // index in the sparsity representation
        while (ind < nnz){
            if (row_ind[ind] == row_offset + nz_rows_[j] && 
                    col_ind[ind] == col_offset + nz_cols_[j]){
                new_nz = false;
                break;
            }
            ind++;
        }

        // if it is a new nonzero, append the sparsity representation
        if (new_nz){
            row_ind[nnz] = row_offset + nz_rows_[j];
            col_ind[nnz] = col_offset + nz_cols_[j];
            if (hessian){
                bookkeeper_.hess_nnz_++;
                nnz++;
            } else {
                bookkeeper_.jac_nnz_++;
                nnz++;
            }
        }

        // note where this nonzero has been added in the sparsity so that
        // the NLPInterface object can easily add numerical values
        nz_indices_[j] = ind;
    }
    return nnz_added;
};

void AdaptiveNLP::addCompleteJacobianContribution(const Sparsity&sp,
                                        std::vector<int> k_indices,
                                        bool includes_time,
                                        int row_offset,
                                        std::vector<std::optional<int>>& 
                                                            bookkeeper_vector){
    // NOTE: this function will always automatically add nonzeros in the
    // ordering imposed by NLPInterface::getNonzeroIndices so no special care
    // needs to be taken here
    if (sp.nnz() == 0){ return;};    

    // sparsity submatarix
    Sparsity sp_sub;

    // rows to be selected to get the submatrix
    std::vector<casadi_int> selected_rows_jac(sp.size1());
    for (int i = 0; i < sp.size1(); i++){ selected_rows_jac[i] = i;}
    
    // pointer to the index of the bookkeeper vector that we can write to
    int nz_ind_ptr = 0;

    // if this sparsity also includes a time-contribution, treat the first
    // column of this sparsity first
    if (includes_time){
        std::vector<casadi_int> time_col(1, 0);
        sp_sub = sp.sub(selected_rows_jac, time_col, mapping_);
        addElementaryContribution(sp_sub, row_offset, 0, false,
                                  bookkeeper_.jac_row_ind_,
                                  bookkeeper_.jac_col_ind_);
        std::copy(nz_indices_.begin(), nz_indices_.begin() + sp_sub.nnz(), 
                  bookkeeper_vector.begin() + nz_ind_ptr);
        nz_ind_ptr += sp_sub.nnz();
    }

    // process the other columns of the sparsity
    int col_offset;
    for (int i = 0; i < k_indices.size(); i++){
        col_offset = (nx_+nu_)*k_indices[i] - 
                     nu_*(k_indices[i] > bookkeeper_.final_ind_) + 
                     free_time_;

        // get the subsparsity corresponding to the correct k_indices
        if (k_indices[i] == bookkeeper_.final_ind_ || 
                sp.size2() == nx_ + includes_time ||
                (k_indices.size() > 1 && i == k_indices.size()-1)){
            col_ptr_ = &selected_cols_nx_;
        } else { col_ptr_ = &selected_cols_nxnu_;}
        
        for (int j = 0; j < (*col_ptr_).size(); j++){
            (*col_ptr_)[j] = (nx_+nu_)*i + includes_time + j;
        }
        sp_sub = sp.sub(selected_rows_jac, *col_ptr_, mapping_);

        if (sp_sub.nnz() > 0){
            // process the subsparsity
            addElementaryContribution(sp_sub, row_offset, col_offset, false,
                                      bookkeeper_.jac_row_ind_,
                                      bookkeeper_.jac_col_ind_);

            // update bookkeeper information
            std::copy(nz_indices_.begin(), nz_indices_.begin() + sp_sub.nnz(), 
                      bookkeeper_vector.begin() + nz_ind_ptr);
            nz_ind_ptr += sp_sub.nnz();
        }
    }
};

void AdaptiveNLP::addCompleteHessianContribution(const Sparsity& sp,
                                    std::vector<int> k_indices,
                                    bool includes_time,
                                    std::vector<std::optional<int>>& 
                                                            bookkeeper_vector){
    // Note on implementation:
    //      First, the first column is processed if this sparsity contains a
    //      time-dependency
    //      Next, the main body of the sparsity is considered.
    //          - The sparsity is split in elementary sparsities (using 
    //            k_indices). The lower triangular elementary sparsities are
    //            then added separately.
    //          - For every elementary sparsity, the sparsity indices of the
    //            nonzero elements are stored in the addElementaryContribution
    //            function. However, the bookkeeper_vector expects to find the
    //            sparsity indices of the nonzero elements in a column-wise
    //            order. So for every elementary sparsity, the sparsity indices
    //            of the nonzero elements have to be written at specific
    //            locations of the bookkeeper_vector.
    //            The vectors 'nnz_per_col_written' and 
    //            'nnz_per_total_cumulative' are defined to help with this.

    if (sp.nnz() == 0){ return;};

    // sparsity submatarix
    Sparsity sp_sub;

    // rows to be selected
    for (int i = 0; i < nx_ + nu_; i++){ selected_rows_nxnu_[i] = i;}

    // keep track of how many nonzeros are written in every column
    std::vector<int> nnz_per_col_written(sp.size2(), 0);

    // For every column, compute how many (lower triangular) nonzeros there
    // are in the sparsity in the columns before
    std::vector<int> nnz_per_col_total_cumulative(sp.size2(), 0);
    int nnz_added = getNonzeroIndices(sp, true);
    for (int i = 0; i < nnz_added; i++){
        for (int j = nz_cols_[i]+1; j < sp.size2(); j++){
            nnz_per_col_total_cumulative[j]++;
        }
    }

    // colind of subsparsity
    std::vector<casadi_int> colind_sub;

    // offsets of the subsparsities to add in the complete hessian sparsity
    int row_offset;
    int col_offset = 0;

    // If the sparsity includes a time dependency, process the first column
    // of the hessian first (related to time)
    if (includes_time){
        // top-left element
        std::vector<casadi_int> time_col(1, 0);
        std::vector<casadi_int> time_row(1, 0);
        sp_sub = sp.sub(time_row, time_col, mapping_);
        if (sp_sub.has_nz(0,0)){
            addElementaryContribution(sp_sub, 0, 0, true,
                                      bookkeeper_.hess_row_ind_,
                                      bookkeeper_.hess_col_ind_);
            bookkeeper_vector[0] = nz_indices_[0];
            nnz_per_col_written[0]++;
        }
        
        // process the first column for dfifferent time steps
        for (int k = 0; k < k_indices.size(); k++){
            // check which rows to select
            if (k_indices[k] == bookkeeper_.final_ind_ || 
                    sp.size2() == nx_+includes_time ||
                    (k_indices.size() > 1 && k == k_indices.size()-1)){
                row_ptr_ = &selected_rows_nx_;
            } else { row_ptr_ = &selected_rows_nxnu_;}
            
            // initialize the row selection
            for (int i = 0; i < (*row_ptr_).size(); i++){
               (*row_ptr_)[i] = 1 + (nx_+nu_)*k + i;
            }

            // get the subsparsity
            sp_sub = sp.sub(*row_ptr_, time_col, mapping_);
            
            // process the subsparsity
            if (sp_sub.nnz() > 0){
                row_offset = 1 + (nx_+nu_)*k_indices[k] - 
                            (k_indices[k] > bookkeeper_.final_ind_)*nu_;
                addElementaryContribution(sp_sub, row_offset, col_offset, true,
                                          bookkeeper_.hess_row_ind_, 
                                          bookkeeper_.hess_col_ind_);
                
                // update the bookkeeping information (keep track of the
                // sparsity indices of the nonzero elements)
                for (int j = 0; j < sp_sub.nnz(); j++){
                    bookkeeper_vector[nnz_per_col_written[0]] = nz_indices_[j];
                    nnz_per_col_written[0]++;
                }
            }
        }
    }

    // add lower triangular blocks
    bool use_offset;
    for (int k = 0; k < k_indices.size(); k++){
        // select correct (numnber of) columns
        if (k_indices[k] == bookkeeper_.final_ind_ || 
                sp.size2() == nx_+includes_time ||
                (k_indices.size() > 1 && k == k_indices.size()-1)){
            col_ptr_ = &selected_cols_nx_;
        } else { col_ptr_ = &selected_cols_nxnu_;}

        for (int i = 0; i < (*col_ptr_).size(); i++){
            (*col_ptr_)[i] = includes_time + (nx_+nu_)*k + i;
        }

        for (int j = k; j < k_indices.size(); j++){
            // select correct (number of) rows and get the subsparsity
            if (k_indices[j] == bookkeeper_.final_ind_ || 
                    sp.size2() == nx_+includes_time ||
                    (k_indices.size() > 1 && j == k_indices.size()-1)){
                row_ptr_ = &selected_rows_nx_;
            } else { row_ptr_ = &selected_rows_nxnu_;}

            for (int i = 0; i < (*row_ptr_).size(); i++){
               (*row_ptr_)[i] = includes_time + (nx_+nu_)*j + i;
            }

            sp_sub = sp.sub(*row_ptr_, *col_ptr_, mapping_);

            // process the subsparsity
            if (sp_sub.nnz() != 0){
                // locate subsparsity in complete hessian matrix
                row_offset = (nx_ + nu_)*k_indices[j] - 
                             (k_indices[j] > bookkeeper_.final_ind_)*nu_ + 
                              free_time_;
                col_offset = (nx_ + nu_)*k_indices[k] - 
                             (k_indices[k] > bookkeeper_.final_ind_)*nu_ + 
                              free_time_;

                // add sparsity to hessian matrix
                nnz_added = addElementaryContribution(sp_sub, row_offset,
                                                    col_offset, true,
                                                    bookkeeper_.hess_row_ind_,
                                                    bookkeeper_.hess_col_ind_);

                // update the bookkeeping information (keep track of the
                // sparsity indices of the nonzero elements)
                for (int i = 0; i < nnz_added; i++){
                    // make sure to write the sparsity index of the nonzero
                    // in the correct place in the bookkeeper_vector
                    bookkeeper_vector[
                        nnz_per_col_total_cumulative[(*col_ptr_)[0] + 
                                                     nz_cols_[i]] +
                        nnz_per_col_written[(*col_ptr_)[0] + nz_cols_[i]]] = 
                            nz_indices_[i];
                    nnz_per_col_written[(*col_ptr_)[0] + nz_cols_[i]]++;
                }
            } 
        }
    }
}

bool AdaptiveNLP::getFirstConstraintAndInstanceInd(int k, int& constraint_ind,
                                                   int& instance_ind){
    for (int i = 0; i < blocks_->nb_extra_blocks_; i++){
        for (int j = 0; j < bookkeeper_.max_nb_extra_instances_; j++){
            if (bookkeeper_.ind_g_e_.count({k, i, j}) > 0){
                constraint_ind = i;
                instance_ind = j;
                return true;
            }
        }
    }
    // no extra constraint found
    return false;
};

bool AdaptiveNLP::getFirstInstanceInd(int k, int constraint_ind, 
                                      int& instance_ind){
    for (int j = 0; j < bookkeeper_.max_nb_extra_instances_; j++){
        if (bookkeeper_.ind_g_e_.count({k, constraint_ind, j}) > 0){
            instance_ind = j;
            return true;
        }
    }
    // no instance found
    return false;
};

bool AdaptiveNLP::getFreeInstanceSlot(int k, int constraint_ind, 
                                      int& instance_ind){
    for (int i = 0; i < bookkeeper_.max_nb_extra_instances_; i++){
        if (bookkeeper_.ind_g_e_.count({k, constraint_ind, i}) == 0){
            instance_ind = i;
            return true;
        }
    }
    return false;
};