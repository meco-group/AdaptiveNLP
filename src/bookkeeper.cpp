#include "../include/bookkeeper.hpp"
#include <cassert>
#include <iostream>

using namespace std;

// constructor
Bookkeeper::Bookkeeper(BuildingBlocks& blocks, int Nmax, 
                       int max_nb_constraints,
                       int max_nb_extra_instances)
                       :Nmax_(Nmax), nx_(blocks.nx_), nu_(blocks.nu_),
                       max_nb_constraints_(max_nb_constraints),
                       nb_extra_blocks_(blocks.nb_extra_blocks_),
                       nb_disc_blocks_(blocks.nb_disc_blocks_),
                       max_nb_extra_instances_(max_nb_extra_instances),
                       free_time_(blocks.free_time_){
    // Initialize lower and upper bound vectors
    lb_container_ = std::vector<double>(max_nb_constraints_);
    ub_container_ = std::vector<double>(max_nb_constraints_);

    // Initialize number of extra constraints applied;
    nb_g_extra_applied_ = std::vector<int>(nb_extra_blocks_);

    // Initialize jacobian sparsity
    jac_row_ind_ = std::vector<std::optional<int>>
            (max_nb_constraints_*(nx_*(Nmax_+1) + nu_*Nmax_));
    jac_col_ind_ = std::vector<std::optional<int>>
            (max_nb_constraints_*(nx_*(Nmax_+1) + nu_*Nmax_));
    jac_nnz_ = 0;

    // Initialize hessian sparsity
    hess_row_ind_ = std::vector<std::optional<int>>
            ((nx_*(Nmax_+1) + nu_*Nmax_)*(nx_*(Nmax_+1) + nu_*Nmax_)/2);
    hess_col_ind_ = std::vector<std::optional<int>>
            ((nx_*(Nmax_+1) + nu_*Nmax_)*(nx_*(Nmax_+1) + nu_*Nmax_)/2);
    hess_nnz_ = 0;

    // Initialize indices in the sparsity of jacobian nonzero contributions
    jac_nz_ind_g0_ = std::vector<std::optional<int>>(blocks.get_nb_jac_nnz(0));
    jac_nz_ind_gT_ = std::vector<std::optional<int>>(blocks.get_nb_jac_nnz(1));
    jac_nz_ind_g_fixed_ = std::vector<std::vector<std::optional<int>>>(Nmax_+1, 
        std::vector<std::optional<int>>(blocks.get_nb_jac_nnz(2)));
    jac_nz_ind_g_disc_ = std::vector<std::vector<std::optional<int>>>
        (Nmax_+1, std::vector<std::optional<int>>(blocks.max_disc_J_nnz_));
    jac_nz_ind_g_extra_ = 
        std::vector<std::vector<std::vector<std::vector<std::optional<int>>>>>
            (Nmax_+1, std::vector<std::vector<std::vector<std::optional<int>>>>
                (nb_extra_blocks_, std::vector<std::vector<std::optional<int>>>
                    (max_nb_extra_instances_)));
    for (int k = 0; k < Nmax_+1; k++){
        for (int j = 0; j < nb_extra_blocks_; j++){
            for (int l = 0; l < max_nb_extra_instances_; l++){
                jac_nz_ind_g_extra_[k][j][l] = 
                    std::vector<std::optional<int>>(blocks.
                        get_nb_jac_nnz_extra(j));
            }
        }
    }

    // Initialize indices in the sparsity of hessian nonzero contributions
    hess_nz_ind_Phi_0_ = std::vector<std::optional<int>>
                            (blocks.get_nb_hess_nnz(0));
    hess_nz_ind_phi_ = std::vector<std::vector<std::optional<int>>>(Nmax_+1, 
        std::vector<std::optional<int>>(blocks.get_nb_hess_nnz(1)));
    hess_nz_ind_Phi_f_ = std::vector<std::optional<int>>
                            (blocks.get_nb_hess_nnz(2));
    hess_nz_ind_g0_ = std::vector<std::optional<int>>
                            (blocks.get_nb_hess_nnz(3));
    hess_nz_ind_gT_ = std::vector<std::optional<int>>
                            (blocks.get_nb_hess_nnz(4));
    hess_nz_ind_g_fixed_ = std::vector<std::vector<std::optional<int>>>
        (Nmax_+1, std::vector<std::optional<int>>(blocks.get_nb_hess_nnz(5)));
    hess_nz_ind_g_disc_ = std::vector<std::vector<std::optional<int>>>
        (Nmax+1, std::vector<std::optional<int>>(blocks.max_disc_H_nnz_));
    hess_nz_ind_g_extra_ = 
        std::vector<std::vector<std::vector<std::vector<std::optional<int>>>>>
        (Nmax_+1, std::vector<std::vector<std::vector<std::optional<int>>>>
                (nb_extra_blocks_, std::vector<std::vector<std::optional<int>>>
                    (max_nb_extra_instances_)));
    for (int k = 0; k < Nmax_+1; k++){
        for (int j = 0; j < nb_extra_blocks_; j++){
            for (int l = 0; l < max_nb_extra_instances_; l++){
                hess_nz_ind_g_extra_[k][j][l] = 
                    std::vector<std::optional<int>>(blocks.
                        get_nb_hess_nnz_extra(j));
            }
        }
    }

    // Initialize constraint indices
    ind_g_f_ = std::vector<std::optional<int>>(Nmax_+1);
    ind_g_d_ = std::vector<std::optional<int>>(Nmax_+1);
    nb_g_d_ = std::vector<std::optional<int>>(Nmax_+1);
    ind_disc_block_ = std::vector<std::optional<int>>(Nmax_+1);
    ind_x_d_ = std::vector<std::optional<int>>(Nmax_+1);

    next_free_ind_ = 0;
    final_ind_ = 0;
    next_ind_ = std::vector<std::optional<int>>(Nmax_+1);
    time_from_ind_ = std::vector<std::optional<double>>(Nmax_+1);

    std::cout<<"bookkeeper is constructed"<<std::endl;
};

double Bookkeeper::dt(int k, double T){
    assert (time_from_ind_[k].has_value());
    assert (next_ind_[k].has_value());
    assert (time_from_ind_[next_ind_[k].value()].has_value());
    if (free_time_){
        return T*(time_from_ind_.at(next_ind_.at(k).value()).value() - \
                  time_from_ind_.at(k).value()); 
    } else {
        return (time_from_ind_.at(next_ind_.at(k).value()).value() - \
                  time_from_ind_.at(k).value());
    }
}

int Bookkeeper::next(int k, int k_break, int k_attach){
    if (k == k_break){ return k_attach;} 
    else { return next(k);}
}

int Bookkeeper::next(int k){
    assert (next_ind_.at(k).has_value());
    return next_ind_.at(k).value();
}

int Bookkeeper::previous(int k){
    int previous_ind = 0;
    while (next(previous_ind) != k){
        previous_ind = next(previous_ind);
    }
    return previous_ind;
}