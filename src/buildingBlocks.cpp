#include "../include/buildingBlocks.hpp"
#include "../include/customExceptions.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <optional>

using namespace casadi;
using namespace std;

// default constructor
BuildingBlocks::BuildingBlocks(){};

// Add problem sizes
void BuildingBlocks::add_sizes(int nx, int nu, bool free_time){
    if (sizes_initialized_){ throw illegalBlocksChangeException();};
    nx_ = nx;
    nu_ = nu;
    free_time_ = free_time;
    sizes_initialized_ = true;
};

// Add stage cost
void BuildingBlocks::add_phi(Function& phi, Function& phi_J,
                             Function& phi_H, int nb_p){
    if (stage_cost_initialized_){ throw illegalBlocksChangeException();};
    phi_ = phi;
    phi_J_ = phi_J;
    phi_H_ = phi_H;
    nb_p_phi_ = nb_p;
    int temp = std::max({phi.nnz_out(0), phi_J.nnz_out(0), phi_H.nnz_out(0)});
    if (temp > max_nnz_){
        max_nnz_ = temp;
        filter_default_ = std::vector<bool>(max_nnz_, true);
    }

    mapping_phi_J_ = std::vector<int>(phi_J.nnz_out(0));
    initMappingVector(phi_J.sparsity_out(0), mapping_phi_J_);

    filter_phi_H_ = std::vector<bool>(phi_H.nnz_out(0));
    initFilterVector(phi_H.sparsity_out(0), filter_phi_H_);

    phi_J_has_nz_ = phi_J.nnz_out(0) > 0;
    phi_H_has_nz_ = phi_H.nnz_out(0) > 0;

    stage_cost_initialized_ = true;
};

// Add cost to initial state
void BuildingBlocks::add_Phi_0(Function& Phi_0, Function& Phi_0_J,
                               Function& Phi_0_H, int nb_p){
    if (initial_cost_initialized_){ throw illegalBlocksChangeException();};
    Phi_0_ = Phi_0;
    Phi_0_J_ = Phi_0_J;
    Phi_0_H_ = Phi_0_H;
    nb_p_Phi_0_ = nb_p;
    int temp = std::max({Phi_0_.nnz_out(0), Phi_0_J_.nnz_out(0), 
                        Phi_0_H_.nnz_out(0)});
    if (temp > max_nnz_){
        max_nnz_ = temp;
        filter_default_ = std::vector<bool>(max_nnz_, true);
    }

    mapping_Phi_0_J_ = std::vector<int>(Phi_0_J_.nnz_out(0));
    initMappingVector(Phi_0_J_.sparsity_out(0), mapping_Phi_0_J_);

    filter_Phi_0_H_ = std::vector<bool>(Phi_0_H_.nnz_out(0));
    initFilterVector(Phi_0_H_.sparsity_out(0), filter_Phi_0_H_);

    Phi_0_J_has_nz_ = Phi_0_J_.nnz_out(0) > 0;
    Phi_0_H_has_nz_ = Phi_0_H_.nnz_out(0) > 0;

    initial_cost_initialized_ = true;
};

// Add cost to final state
void BuildingBlocks::add_Phi_f(Function& Phi_f, Function& Phi_f_J,
                               Function& Phi_f_H, int nb_p){
    if (final_cost_initialized_){ throw illegalBlocksChangeException();};
    Phi_f_ = Phi_f;
    Phi_f_J_ = Phi_f_J;
    Phi_f_H_ = Phi_f_H;
    nb_p_Phi_f_ = nb_p;
    int temp = std::max({Phi_f_.nnz_out(0), Phi_f_J_.nnz_out(0), 
                        Phi_f_H_.nnz_out(0)});
    if (temp > max_nnz_){
        max_nnz_ = temp;
        filter_default_ = std::vector<bool>(max_nnz_, true);
    }

    mapping_Phi_f_J_ = std::vector<int>(Phi_f_J_.nnz_out(0));
    initMappingVector(Phi_f_J_.sparsity_out(0), mapping_Phi_f_J_);

    filter_Phi_f_H_ = std::vector<bool>(Phi_f_H_.nnz_out(0));
    initFilterVector(Phi_f_H_.sparsity_out(0), filter_Phi_f_H_);

    Phi_f_J_has_nz_ = Phi_f_J_.nnz_out(0) > 0;
    Phi_f_H_has_nz_ = Phi_f_H_.nnz_out(0) > 0;

    final_cost_initialized_ = true;
};

// Add initial conditions
void BuildingBlocks::add_g0(Function& g0, Function& g0_J, Function& g0_H,
                            std::vector<double>& lb, std::vector<double>& ub,
                            int nb_p){
    if (initial_constraint_initialized_){
        throw illegalBlocksChangeException();};
    g0_ = g0;
    g0_J_ = g0_J;
    g0_H_ = g0_H;
    nb_g0_ = g0_.size1_out(0);
    if (nb_g0_ > max_nb_g_){
        max_nb_g_ = nb_g0_;
    }
    assert (lb.size() == nb_g0_);
    assert (ub.size() == nb_g0_);
    g0_lb_ = lb;
    g0_ub_ = ub;
    nb_p_g0_ = nb_p;
    int temp = std::max({g0_.nnz_out(0), g0_J_.nnz_out(0), g0_H_.nnz_out(0)});
    if (temp > max_nnz_){
        max_nnz_ = temp;
        filter_default_ = std::vector<bool>(max_nnz_, true);
    }

    filter_g0_H_ = std::vector<bool>(g0_H_.nnz_out(0));
    initFilterVector(g0_H_.sparsity_out(0), filter_g0_H_);

    g0_J_has_nz_ = g0_J_.nnz_out(0) > 0;
    g0_H_has_nz_ = g0_H_.nnz_out(0) > 0;

    initial_constraint_initialized_ = true;
};

// Add final conditions
void BuildingBlocks::add_gT(Function& gT, Function& gT_J, Function& gT_H,
            std::vector<double>& lb, std::vector<double>& ub, int nb_p){
    if (final_constraint_initialized_){ throw illegalBlocksChangeException();};
    gT_ = gT;
    gT_J_ = gT_J;
    gT_H_ = gT_H;
    nb_gT_ = gT_.size1_out(0);
    if (nb_gT_ > max_nb_g_){
        max_nb_g_ = nb_gT_;
    }
    assert (lb.size() == nb_gT_);
    assert (ub.size() == nb_gT_);
    gT_lb_ = lb;
    gT_ub_ = ub;
    nb_p_gT_ = nb_p;
    int temp = std::max({gT_.nnz_out(0), gT_J_.nnz_out(0), gT_H_.nnz_out(0)});
    if (temp > max_nnz_){
        max_nnz_ = temp;
        filter_default_ = std::vector<bool>(max_nnz_, true);
    }

    filter_gT_H_ = std::vector<bool>(gT_H_.nnz_out(0));
    initFilterVector(gT_H_.sparsity_out(0), filter_gT_H_);

    gT_J_has_nz_ = gT_J_.nnz_out(0) > 0;
    gT_H_has_nz_ = gT_H_.nnz_out(0) > 0;

    final_constraint_initialized_ = true;
};

// Add fixed constraints
void BuildingBlocks::add_g_fixed(Function& g_fixed, Function& g_fixed_J,
                                 Function& g_fixed_H,
                                 std::vector<double>& lb,
                                 std::vector<double>& ub, int nb_p){
    if (fixed_constraint_initialized_){ throw illegalBlocksChangeException();};
    g_fixed_ = g_fixed;
    g_fixed_J_ = g_fixed_J;
    g_fixed_H_ = g_fixed_H;
    nb_g_fixed_ = g_fixed_.size1_out(0);
    if (nb_g_fixed_ > max_nb_g_){
        max_nb_g_ = nb_g_fixed_;
    }
    assert (lb.size() == nb_g_fixed_);
    assert (ub.size() == nb_g_fixed_);
    g_fixed_lb_ = lb;
    g_fixed_ub_ = ub;
    nb_p_g_fixed_ = nb_p;
    int temp = std::max({g_fixed_.nnz_out(0), g_fixed_J_.nnz_out(0),
                        g_fixed_H_.nnz_out(0)});
    if (temp > max_nnz_){
        max_nnz_ = temp;
        filter_default_ = std::vector<bool>(max_nnz_, true);
    }

    filter_g_fixed_H_ = std::vector<bool>(g_fixed_H_.nnz_out(0));
    initFilterVector(g_fixed_H_.sparsity_out(0), filter_g_fixed_H_);

    g_fixed_J_has_nz_ = g_fixed_J_.nnz_out(0) > 0;
    g_fixed_H_has_nz_ = g_fixed_H_.nnz_out(0) > 0;

    fixed_constraint_initialized_ = true;
};

// Add extra constraints
void BuildingBlocks::add_g_extra(std::vector<Function> g_extra,
                    std::vector<Function> g_extra_J,
                    std::vector<Function> g_extra_H,
                    std::vector<std::vector<double>> lb,
                    std::vector<std::vector<double>> ub){
    if (extra_constraint_initialized_){ throw illegalBlocksChangeException();};
    assert (g_extra.size() == g_extra_J.size());
    assert (g_extra.size() == g_extra_H.size());
    g_extra_ = g_extra;
    g_extra_J_ = g_extra_J;
    g_extra_H_ = g_extra_H;
    nb_extra_blocks_ = g_extra_.size();
    nb_g_extra_ = std::vector<int>(nb_extra_blocks_);
    g_extra_J_has_nz_total_ = false;
    g_extra_H_has_nz_total_ = false;
    g_extra_J_has_nz_ = std::vector<bool>(nb_extra_blocks_);
    g_extra_H_has_nz_ = std::vector<bool>(nb_extra_blocks_);
    int temp;
    filter_g_extra_H_ = std::vector<std::vector<bool>>(nb_extra_blocks_);
    for (int i = 0; i < nb_extra_blocks_; i++){
        nb_g_extra_[i] = g_extra_[i].size1_out(0);
        nb_g_extra_total_ += nb_g_extra_[i];
        if (nb_g_extra_[i] > max_nb_g_){
           max_nb_g_ = nb_g_extra_[i];
        }
        temp = std::max({g_extra_[i].nnz_out(0), g_extra_J_[i].nnz_out(0),
                        g_extra_H_[i].nnz_out(0)});

        filter_g_extra_H_[i] = std::vector<bool>(g_extra_H_[i].nnz_out(0));
        initFilterVector(g_extra_H_[i].sparsity_out(0), filter_g_extra_H_[i]);
        
        g_extra_J_has_nz_[i] = g_extra_J_[i].nnz_out(0) > 0;
        g_extra_H_has_nz_[i] = g_extra_H_[i].nnz_out(0) > 0;
        g_extra_J_has_nz_total_ = g_extra_J_has_nz_total_ || 
                                  g_extra_J_has_nz_[i];
        g_extra_H_has_nz_total_ = g_extra_H_has_nz_total_ || 
                                  g_extra_H_has_nz_[i];

        if (temp > max_nnz_){
            max_nnz_ = temp;
            filter_default_ = std::vector<bool>(max_nnz_, true);
        } 
    }
    g_extra_lb_ = lb;
    g_extra_ub_ = ub;
    nb_p_g_extra_ = std::vector<int>(nb_extra_blocks_, 0);
    extra_constraint_initialized_ = true;
};
                
// Add discretization constraints
void BuildingBlocks::add_g_disc(std::vector<Function> g_disc, 
                                std::vector<Function> g_disc_J,
                                std::vector<Function> g_disc_H){
    if (discretization_constraint_initialized_){ 
        throw illegalBlocksChangeException();};
    assert (g_disc.size() == g_disc_J.size());
    assert (g_disc.size() == g_disc_H.size());
    g_disc_ = g_disc;
    g_disc_J_ = g_disc_J;
    g_disc_H_ = g_disc_H;
    nb_disc_blocks_ = g_disc_.size();
    nb_g_disc_ = std::vector<int>(g_disc_.size());
    g_disc_J_has_nz_total_ = false;
    g_disc_H_has_nz_total_ = false;
    g_disc_J_has_nz_ = std::vector<bool>(g_disc_.size());
    g_disc_H_has_nz_ = std::vector<bool>(g_disc_.size());
    max_disc_J_nnz_ = 0;
    max_disc_H_nnz_ = 0;
    int temp;
    filter_g_disc_H_ = std::vector<std::vector<bool>>(g_disc_.size());
    for (int i = 0; i < g_disc_.size(); i++){
        nb_g_disc_[i] = g_disc_[i].size1_out(0);
        if (nb_g_disc_[i] > max_nb_g_){
            max_nb_g_ = nb_g_disc_[i];
        }
        max_disc_J_nnz_ = std::max(max_disc_J_nnz_, 
                                   int(g_disc_J_[i].nnz_out(0)));
        max_disc_H_nnz_ = std::max(max_disc_H_nnz_, 
                                   int(g_disc_H_[i].nnz_out(0)));

        if (g_disc_[i].sparsity_in(0).size2() > max_nb_steps_){
            max_nb_steps_ = g_disc_[i].sparsity_in(0).size2();
        }
        int temp = std::max({g_disc_[i].nnz_out(0), g_disc_J_[i].nnz_out(0),
                            g_disc_H_[i].nnz_out(0)});

        filter_g_disc_H_[i] = std::vector<bool>(g_disc_H_[i].nnz_out(0));
        initFilterVector(g_disc_H_[i].sparsity_out(0), filter_g_disc_H_[i]);

        g_disc_J_has_nz_[i] = g_disc_J_[i].nnz_out(0) > 0;
        g_disc_H_has_nz_[i] = g_disc_H_[i].nnz_out(0) > 0;
        g_disc_J_has_nz_total_ = g_disc_J_has_nz_total_ || 
                                 g_disc_J_has_nz_[i];
        g_disc_H_has_nz_total_ = g_disc_H_has_nz_total_ || 
                                 g_disc_H_has_nz_[i];
        
        if (temp > max_nnz_){
            max_nnz_ = temp;
            filter_default_ = std::vector<bool>(max_nnz_, true);
        }
    }

    nb_steps_mapping_ = std::vector<std::optional<int>>(max_nb_steps_+1);
    for (int i = 0; i < g_disc_.size(); i++){
        if (nb_steps_mapping_[g_disc_[i].size2_in(0)].has_value()){
            throw illegalBlocksChangeException();
        }
        nb_steps_mapping_[g_disc_[i].size2_in(0)] = i;
    }

    nb_p_g_disc_ = std::vector<int>(nb_disc_blocks_, 0);
    discretization_constraint_initialized_ = true;
};

// Add initial guesses
void BuildingBlocks::add_inits(Function& x_init, Function& u_init){
    if (init_functions_initialized_){ throw illegalBlocksChangeException();};
    x_init_ = x_init;
    u_init_ = u_init;
    t_init_ = -1;
    init_functions_initialized_ = true;
};
void BuildingBlocks::add_inits(Function& x_init, Function& u_init,
                               double t_init){
    if (init_functions_initialized_){ throw illegalBlocksChangeException();};
    x_init_ = x_init;
    u_init_ = u_init;
    t_init_ = t_init;
    init_functions_initialized_ = true;
};

// check complete collection of building blocks
bool BuildingBlocks::check_complete(){
    return sizes_initialized_ && stage_cost_initialized_ && 
        initial_cost_initialized_ && final_cost_initialized_ && 
        initial_constraint_initialized_ && final_constraint_initialized_ && 
        fixed_constraint_initialized_ && extra_constraint_initialized_ &&
        discretization_constraint_initialized_ && init_functions_initialized_;
};

bool BuildingBlocks::check_problem_sizes(int nx_check, int nu_check){
    return nx_ == nx_check && nu_ == nu_check;
};

int BuildingBlocks::get_nb_jac_nnz(int block_nb){
    switch(block_nb){
        case 0:
            return g0_J_.nnz_out(0);
        case 1:
            return gT_J_.nnz_out(0);
        case 2:
            return g_fixed_J_.nnz_out(0);
        default:
            return -1;
    }
}

int BuildingBlocks::get_nb_jac_nnz_extra(int extra_nb){
    return g_extra_J_[extra_nb].nnz_out(0);
}

int BuildingBlocks::get_nb_jac_nnz_disc(int disc_nb){
    return g_disc_J_[disc_nb].nnz_out(0);
}

int BuildingBlocks::get_nb_hess_nnz(int block_nb){
    switch(block_nb){
        case 0:
            return Phi_0_H_.nnz_out(0);
        case 1:
            return phi_H_.nnz_out(0);
        case 2:
            return Phi_f_H_.nnz_out(0);
        case 3:
            return g0_H_.nnz_out(0);
        case 4:
            return gT_H_.nnz_out(0);
        case 5:
            return g_fixed_H_.nnz_out(0);
        default:
            return -1;
    }
}

int BuildingBlocks::get_nb_hess_nnz_extra(int extra_nb){
    return g_extra_H_[extra_nb].nnz_out(0);
}

int BuildingBlocks::get_nb_hess_nnz_disc(int disc_nb){
    return g_disc_H_[disc_nb].nnz_out(0);
}

void BuildingBlocks::initMappingVector(const Sparsity& sp,
                                       std::vector<int>& mapping){
    std::vector<casadi_int> colind = sp.get_colind();

	int idx_ptr = 0;
	for (int col_ptr = 1; col_ptr < sp.size2()+1; col_ptr++){
        if (colind[col_ptr] > colind[col_ptr-1]){
            mapping[idx_ptr] = col_ptr-1;
            idx_ptr++;
        }
	}
};

void BuildingBlocks::initFilterVector(const Sparsity& sp,
                                      std::vector<bool>& filter){
    std::vector<casadi_int> colind = sp.get_colind();
	std::vector<casadi_int> rowind = sp.get_row();
	int idx_ptr = 0;
	for (int col_ptr = 1; col_ptr < colind.size(); col_ptr++){
		for (int row_ptr = colind[col_ptr-1]; row_ptr < colind[col_ptr]; 
				row_ptr++){
            if (rowind[row_ptr] >= col_ptr-1){
                filter[idx_ptr] = true;
            } else {
                filter[idx_ptr] = false;
            }
            idx_ptr++;
		}
	}
};