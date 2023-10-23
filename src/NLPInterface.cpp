#include "../include/NLPInterface.hpp"
#include "../include/customExceptions.hpp"
#include <limits>
#include <cmath>
#include <chrono>

using namespace std;
using namespace Ipopt;
using namespace std::chrono;

void printOptionalVectorHere(std::vector<std::optional<int>>& v, int max_ind, bool endline=true){
	for (int i = 0; i < max_ind; i++){
		if (v[i].has_value()){
			std::cout<<std::setw(3)<<std::setfill('0')<<v[i].value()<<" ";
		} else {
			std::cout<<" X  ";
		}
	}
	if (endline){std::cout<<endl;}
}


NLPInterface::NLPInterface(int N, Bookkeeper& bookkeeper,
                    BuildingBlocks& buildingBlocks_, bool diagonalize_)
                    : print_(false), nx_(bookkeeper.nx_),
					  nu_(bookkeeper.nu_), free_time_(bookkeeper.free_time_),
                      Nmax_(bookkeeper.Nmax_), bookkeeper_(&bookkeeper),
					  blocks_(&buildingBlocks_),
					  function_out_(std::max({blocks_->max_nb_g_, nx_+nu_,
											 blocks_->max_nnz_})),
					  res1_(1, &function_out_[0]){
    assert (buildingBlocks_.check_complete());
    assert (buildingBlocks_.check_problem_sizes(nx_, nu_));

    N_ = N;
	init_guess_ = std::vector<double>(nx_*(Nmax_+1) + nu_*Nmax_ + free_time_);
    diagonalize_ = diagonalize_;

	// allocate scratch space
	k_vals_ = std::vector<int>(blocks_->max_nb_steps_);
	dts_ = std::vector<double>(blocks_->max_nb_steps_ - 1);
	xk_ = std::vector<double>(nx_);
	uk_ = std::vector<double>(nu_);
	xx_ = Eigen::MatrixXd(nx_, blocks_->max_nb_steps_);
	uu_ = Eigen::MatrixXd(nu_, blocks_->max_nb_steps_ - 1);

	nz_rows_ = std::vector<int>(blocks_->max_nnz_);
	nz_cols_ = std::vector<int>(blocks_->max_nnz_);

	// initialize clearStructuralZeros scratch space
	x_zeros_ = std::vector<double>(nx_*(Nmax_+1)+nu_*Nmax_ + free_time_, 0.0);
	if (free_time_){x_zeros_[0] = 1.0;}
	lambda_zeros_ = std::vector<double>(bookkeeper_->max_nb_constraints_, 0.0);
	clear_constraints_ = std::vector<double>(bookkeeper_->max_nb_constraints_,
							std::numeric_limits<double>::quiet_NaN());
	clear_jacobian_ = std::vector<double>(bookkeeper_->jac_row_ind_.size(),
							std::numeric_limits<double>::quiet_NaN());
	clear_hessian_ = std::vector<double>(bookkeeper_->hess_row_ind_.size(),
							std::numeric_limits<double>::quiet_NaN());
	g_mapping_ = std::vector<std::optional<int>>(bookkeeper_->
														max_nb_constraints_);
	jac_mapping_ = std::vector<std::optional<int>>(bookkeeper_->
														jac_row_ind_.size());
	hess_mapping_ = std::vector<std::optional<int>>(bookkeeper_->
														hess_row_ind_.size());

	k_in_use_ = std::vector<double>(Nmax_);
	k_mapping_ = std::vector<std::optional<int>>(Nmax_);
	col_in_use_ = std::vector<double>(nx_*(Nmax_+1)+nu_*Nmax_ + free_time_);
	col_mapping_ = std::vector<std::optional<int>>(nx_*(Nmax_+1)+nu_*Nmax_ + 
														free_time_);

	// initialize parameters
	p_Phi_0_ = std::vector<double>(blocks_->nb_p_Phi_0_, 0.0);
    p_Phi_f_ = std::vector<double>(blocks_->nb_p_Phi_f_, 0.0);
    p_phi_ = std::vector<double>(blocks_->nb_p_phi_, 0.0);
    p_g0_ = std::vector<double>(blocks_->nb_p_g0_, 0.0);
    p_gT_ = std::vector<double>(blocks_->nb_p_gT_, 0.0);
    p_g_fixed_ = std::vector<double>(blocks_->nb_p_g_fixed_, 0.0);
    p_g_disc_ = std::vector<std::vector<double>>
		(blocks_->nb_disc_blocks_);
	for (int i = 0; i < blocks_->nb_disc_blocks_; i++){
		p_g_disc_[i] = std::vector<double>(blocks_->nb_p_g_disc_[i], 0.0);
	}
	good_p_g_disc_ = std::vector<bool>(blocks_->max_nb_steps_+1, false);
    p_g_extra_ = std::vector<std::vector<std::vector<std::vector<double>>>>
                	(bookkeeper_->Nmax_+1,
                	std::vector<std::vector<std::vector<double>>>(
                    	blocks_->nb_extra_blocks_,
                		std::vector<std::vector<double>>(
                      		bookkeeper_->max_nb_extra_instances_)));
    for (int k = 0; k < bookkeeper_->Nmax_+1; k++){
    	for (int j = 0; j < blocks_->nb_extra_blocks_; j++){
    		for (int l = 0; l < bookkeeper_->max_nb_extra_instances_; l++){
          		p_g_extra_[k][j][l] = std::vector<double>(blocks_->
															nb_p_g_extra_[j]);
        }
      }
    }

	sol_x_ = std::vector<double>(nx_*(N_+1) + nu_*N_ + free_time_);
	sol_lambda_ = std::vector<double>(bookkeeper_->max_nb_constraints_);

	std::cout<<"interface is constructed"<<std::endl;
};

bool NLPInterface::get_nlp_info(
	Index&          n,
	Index&          m,
	Index&          nnz_jac_g,
	Index&          nnz_h_lag,
	IndexStyleEnum& index_style
){
	if(print_){std::cout<<"entering "<<__func__<<std::endl;}
	try{
	n = nx_*(N_+1) + nu_*N_ + free_time_;
	m = bookkeeper_->next_free_ind_;
	nnz_jac_g = bookkeeper_->jac_nnz_;
	nnz_h_lag = bookkeeper_->hess_nnz_;
	index_style = C_STYLE;

	} catch (const std::exception &exc){std::cerr<<exc.what()<<endl; 
             std::exit(1);}
	if(print_){std::cout<<"exit "<<__func__<<std::endl;}
	return true;
};

bool NLPInterface::get_bounds_info(
	Index   n,
	Number* x_l,
	Number* x_u,
	Index   m,
	Number* g_l,
	Number* g_u
){
	// auto start = high_resolution_clock::now();
	if(print_){std::cout<<"entering "<<__func__<<std::endl;}
	try{
	for (int i = 0; i < nx_*(N_+1) + nu_*N_ + free_time_; i++){
		x_l[i] = -1.0e39;
		x_u[i] =  1.0e39;
	}

	for (int i = 0; i < bookkeeper_->next_free_ind_; i++){
		g_l[i] = bookkeeper_->lb_container_[i];
	}
	for (int i = 0; i < bookkeeper_->next_free_ind_; i++){
		g_u[i] = bookkeeper_->ub_container_[i];
	}

	} catch (const std::exception &exc){std::cerr<<exc.what()<<endl; 
             std::exit(1);}
	if(print_){std::cout<<"exit "<<__func__<<std::endl;}
	// auto stop = high_resolution_clock::now();
	// cout<<"nb of nanoseconds: "<<double(duration_cast<nanoseconds>(stop-start).count())<<endl;
	// nanosecond_counter += double(duration_cast<nanoseconds>(stop-start).count());
	return true;
};

bool NLPInterface::get_starting_point(
	Index   n,
	bool    init_x,
	Number* x,
	bool    init_z,
	Number* z_L,
	Number* z_U,
	Index   m,
	bool    init_lambda,
	Number* lambda
){
	if(print_){std::cout<<"entering "<<__func__<<std::endl;}
	try{
	if (init_x){
		if (initial_guess_available_){
			std::copy(init_guess_.begin(), 
					init_guess_.begin() + nx_*(N_+1) + nu_*N_ + free_time_, x);
		} else {
			int init_pointer = 0;
			if (free_time_){
				x[init_pointer] = blocks_->t_init_;
				init_pointer++;
			}
			for (int k = 0; k <= N_; k++){
				arg1_[0] = &bookkeeper_->time_from_ind_[k].value();
				blocks_->x_init_(arg1_, res1_);
				
				std::copy(function_out_.begin(), function_out_.begin()+nx_,
						  &x[init_pointer]);
				init_pointer += nx_;

				if (k != bookkeeper_->final_ind_){
					blocks_->u_init_(arg1_, res1_);

					std::copy(function_out_.begin(), function_out_.begin()+nu_,
							  &x[init_pointer]);
					init_pointer += nu_;
				}
			}
		}
	}

	// TODO: initialize z and lambda

	} catch (const std::exception &exc){std::cerr<<exc.what()<<endl; 
             std::exit(1);}
	if(print_){std::cout<<"exit "<<__func__<<std::endl;}
	return true;
};

bool NLPInterface::eval_f(
	Index         n,
	const Number* x,
	bool          new_x,
	Number&       obj_value
){
	if(print_){std::cout<<"entering "<<__func__<<std::endl;}
	try{

	if (free_time_){	
		arg3_[0] = &xk_[0];
		arg3_[1] = &x[0];
		
		// add initial cost
		get_x(x, 0);
		arg3_[2] = &p_Phi_0_[0];
		blocks_->Phi_0_(arg3_, res1_);
		obj_value = function_out_[0];
		
		// add final cost
		get_x(x, bookkeeper_->final_ind_);
		arg3_[2] = &p_Phi_0_[0];
		blocks_->Phi_f_(arg3_, res1_);
		obj_value += function_out_[0];

		// add stage cost
		double dt;
		arg5_[0] = &xk_[0]; arg5_[1] = &uk_[0]; arg5_[2] = &dt;
		arg5_[3] = &x[0]; arg5_[4] = &p_phi_[0];
		for (int k = 0; k < N_+1; k++){
			if (k != bookkeeper_->final_ind_){
				get_x(x, k);
				get_u(x, k);
				dt = bookkeeper_->dt(k, x[0]);
				blocks_->phi_(arg5_, res1_);
				obj_value += function_out_[0];
			}
		}
	} else {
		arg2_[0] = &xk_[0];
		
		// add initial cost
		get_x(x, 0);
		arg2_[1] = &p_Phi_0_[0];
		blocks_->Phi_0_(arg2_, res1_);
		obj_value = function_out_[0];

		// add final cost
		get_x(x, bookkeeper_->final_ind_);
		arg2_[1] = &p_Phi_f_[0];
		blocks_->Phi_f_(arg2_, res1_);
		obj_value += function_out_[0];

		// add stage cost
		double dt;
		arg4_[0] = &xk_[0]; arg4_[1] = &uk_[0]; arg4_[2] = &dt;
		arg4_[3] = &p_phi_[0];
		for (int k = 0; k < N_+1; k++){
			if (k != bookkeeper_->final_ind_){
				get_x(x, k);
				get_u(x, k);
				dt = bookkeeper_->dt(k, x[0]);
				blocks_->phi_(arg4_, res1_);
				obj_value += function_out_[0];
			}
		}
	}

	} catch (const std::exception &exc){std::cerr<<exc.what()<<endl; 
             std::exit(1);}
	if(print_){std::cout<<"exit "<<__func__<<std::endl;}
	return true;
};


bool NLPInterface::eval_grad_f(
	Index         n,
	const Number* x,
	bool          new_x,
	Number*       grad_f
){
	if(print_){std::cout<<"entering "<<__func__<<std::endl;}
	try{
	for (int i = 0; i < n; i++){grad_f[i] = 0;}

	if (free_time_){
		arg3_[0] = &xk_[0]; arg3_[1] = &x[0];

		// initial cost
		if (blocks_->Phi_0_J_has_nz_){
			get_x(x, 0);
			arg3_[2] = &p_Phi_0_[0];
			blocks_->Phi_0_J_(arg3_, res1_);
			copyFunctionOutToTarget(grad_f, 0, blocks_->mapping_Phi_0_J_, 
									true);
		}

		// final cost
		if (blocks_->Phi_f_J_has_nz_){
			get_x(x, bookkeeper_->final_ind_);
			arg3_[2] = &p_Phi_f_[0];
			blocks_->Phi_f_J_(arg3_, res1_);
			copyFunctionOutToTarget(grad_f,
									(nx_ + nu_)*bookkeeper_->final_ind_,
									blocks_->mapping_Phi_f_J_, true);
		}

		// stage cost
		if (blocks_->phi_J_has_nz_){
			int offset;
			double dt;
			arg5_[0] = &xk_[0]; arg5_[1] = &uk_[0]; arg5_[2] = &dt; 
			arg5_[3] = &x[0]; arg5_[4] = &p_phi_[0];
			for (int k = 0; k <= N_; k++){
				if (k != bookkeeper_->final_ind_){
					int offset = 1 ? k > bookkeeper_->final_ind_ : 0;
					get_x(x, k); get_u(x, k);
					dt = bookkeeper_->dt(k, x[0]);
					blocks_->phi_J_(arg5_, res1_);
					copyFunctionOutToTarget(grad_f, (nx_+nu_)*k - offset*nu_,
											blocks_->mapping_phi_J_, true);
				}
			}
		}
	} else {
		arg2_[0] = &xk_[0];

		// initial cost
		if (blocks_->Phi_0_J_has_nz_){
			get_x(x, 0);
			arg2_[1] = &p_Phi_0_[0];
			blocks_->Phi_0_J_(arg2_, res1_);
			copyFunctionOutToTarget(grad_f, 0, blocks_->mapping_Phi_0_J_, 
									false);
		}

		// final cost
		if (blocks_->Phi_f_J_has_nz_){
			get_x(x, bookkeeper_->final_ind_);
			arg2_[1] = &p_Phi_f_[0];
			blocks_->Phi_f_J_(arg2_, res1_);
			copyFunctionOutToTarget(grad_f, 
									(nx_ + nu_)*bookkeeper_->final_ind_,
									blocks_->mapping_Phi_f_J_, false);
		}

		// stage cost
		if (blocks_->phi_J_has_nz_){
			int offset;
			double dt;
			arg4_[0] = &xk_[0]; arg4_[1] = &uk_[0]; arg4_[2] = &dt; 
			arg4_[3] = &p_phi_[0];
			for (int k = 0; k <= N_; k++){
				if (k != bookkeeper_->final_ind_){
					int offset = 1 ? k > bookkeeper_->final_ind_ : 0;
					get_x(x, k); get_u(x, k);
					dt = bookkeeper_->dt(k);
					blocks_->phi_J_(arg4_, res1_);
					copyFunctionOutToTarget(grad_f, (nx_+nu_)*k - offset*nu_,
											blocks_->mapping_phi_J_, false);
				}
			}
		}
	}

	} catch (const std::exception &exc){std::cerr<<exc.what();}
	if(print_){std::cout<<"exit "<<__func__<<std::endl;}	
	return true;
};

bool NLPInterface::eval_g(
	Index         n,
	const Number* x,
	bool          new_x,
	Index         m,
	Number*       g
){
	if(print_){std::cout<<"entering "<<__func__<<std::endl;}
	try{
	if (!clearingStructuralZeros_){
		for (int i = 0; i < m; i++){g[i] = 0;}
	}

	// initial constraint
	if (blocks_->nb_g0_ > 0){
		if (free_time_){
			arg3_[0] = &xk_[0]; arg3_[1] = &x[0]; arg3_[2] = &p_g0_[0];
			get_x(x, 0);
			blocks_->g0_(arg3_, res1_);
			for (int i = 0; i < blocks_->nb_g0_; i++){
				g[i] = function_out_[i];
			}
		} else {
			arg2_[0] = &xk_[0]; arg2_[1] = &p_g0_[0];
			get_x(x, 0);
			blocks_->g0_(arg2_, res1_);
			for (int i = 0; i < blocks_->nb_g0_; i++){
				g[i] = function_out_[i];
			}
		}
	}

	// final constraint
	if (blocks_->nb_gT_ > 0){
		get_x(x, bookkeeper_->final_ind_);
		arg2_[0] = &xk_[0]; arg2_[1] = &p_gT_[0];
		blocks_->gT_(arg2_, res1_);
		for (int i = 0; i < blocks_->nb_gT_; i++){
			g[blocks_->nb_g0_ + i] = function_out_[i];
		}
	}

	if (blocks_->nb_g_fixed_ + blocks_->nb_g_extra_total_ > 0){
		arg3_[0] = &xk_[0]; arg3_[1] = &uk_[0];
		for (int k = 0; k <= N_; k++){
			if (k != bookkeeper_->final_ind_){
				// fixed constraints
				get_x(x, k); get_u(x, k);
				if (blocks_->nb_g_fixed_ > 0){
					arg3_[2] = &p_g_fixed_[0];
					blocks_->g_fixed_(arg3_, res1_);
					for (int i = 0; i < blocks_->nb_g_fixed_; i++){
						g[bookkeeper_->ind_g_f_[k].value()+i] = 
							function_out_[i];
					}
				}
			}
		}
	}

	// extra constraints
	arg3_[0] = &xk_[0]; arg3_[1] = &uk_[0];
	for (const auto& pair : bookkeeper_->ind_g_e_){
		get_x(x, pair.first.k); get_u(x, pair.first.k);
		arg3_[2] = &p_g_extra_[pair.first.k][pair.first.j][pair.first.i][0];
		blocks_->g_extra_[pair.first.j](arg3_, res1_);
		for (int i = 0; i < blocks_->nb_g_extra_[pair.first.j]; i++){
			g[pair.second+i] = function_out_[i];
		}
	}

	// dynamical constraints
	int k = 0;
	int k_init;
	
	if (free_time_){
		arg5_[0] = &xx_(0,0); arg5_[1] = &uu_(0,0); arg5_[2] = &dts_[0];
		arg5_[3] = &x[0];
	} else {
		arg4_[0] = &xx_(0,0); arg4_[1] = &uu_(0,0); arg4_[2] = &dts_[0];
	}
	while (k != bookkeeper_->final_ind_){
		k_init = k;
		if (free_time_){
			arg5_[4] = &p_g_disc_[bookkeeper_->ind_disc_block_[k].value()][0];
		} else {
			arg4_[3] = &p_g_disc_[bookkeeper_->ind_disc_block_[k].value()][0];
		}

		// get variable values
		getCollVars(k_init, bookkeeper_->nb_g_d_[k_init].value(), x);

		// get time-steps
		for (int i = 0; i < bookkeeper_->nb_g_d_[k_init].value()-1; i++){
			dts_[i] = bookkeeper_->dt(k_vals_[i]);
		}

		if (free_time_){
			blocks_->g_disc_[bookkeeper_->ind_disc_block_[k].value()]
				(arg5_, res1_);
		} else {
			blocks_->g_disc_[bookkeeper_->ind_disc_block_[k].value()]
				(arg4_, res1_);
		}

		for (int i = 0; i < blocks_->
				nb_g_disc_[bookkeeper_->ind_disc_block_[k].value()]; i++){
			g[bookkeeper_->ind_g_d_[k_init].value()+i] = function_out_[i];
		}

		k = k_vals_[bookkeeper_->nb_g_d_[k_init].value()-1];
	}

	} catch (const std::exception &exc){std::cerr<<exc.what()<<endl; 
             std::exit(1);}
	if(print_){std::cout<<"exit "<<__func__<<std::endl;}
	return true;
};

bool NLPInterface::eval_jac_g(
	Index         n,
	const Number* x,
	bool          new_x,
	Index         m,
	Index         nele_jac,
	Index*        iRow,
	Index*        jCol,
	Number*       values
){
	if(print_){std::cout<<"entering "<<__func__<<std::endl;}
	try{
	if (values == NULL){
		for (int i = 0; i < nele_jac; i++){
			iRow[i] = bookkeeper_->jac_row_ind_[i].value();
			jCol[i] = bookkeeper_->jac_col_ind_[i].value();
		}
	} else {
		if (!clearingStructuralZeros_){
			for (int i = 0; i < nele_jac; i++){values[i] = 0;}
		}

		// initial constraint
		if (blocks_->g0_J_has_nz_){
			if (free_time_){
				arg3_[0] = &xk_[0]; arg3_[1] = &x[0]; arg3_[2] = &p_g0_[0];
				get_x(x, 0);
				blocks_->g0_J_(arg3_, res1_);
				addContribution(values, blocks_->g0_J_.sparsity_out(0).nnz(),
								bookkeeper_->jac_nz_ind_g0_);
			} else {
				arg3_[0] = &xk_[0]; arg2_[1] = &p_g0_[0];
				get_x(x, 0);
				blocks_->g0_J_(arg2_, res1_);
				addContribution(values, blocks_->g0_J_.sparsity_out(0).nnz(),
								bookkeeper_->jac_nz_ind_g0_);
			}
		}

		// final constraint
		if (blocks_->gT_J_has_nz_){
			arg2_[0] = &xk_[0]; arg2_[1] = &p_gT_[0];
			get_x(x, bookkeeper_->final_ind_);
			blocks_->gT_J_(arg2_, res1_);
			addContribution(values, blocks_->gT_J_.sparsity_out(0).nnz(),
							bookkeeper_->jac_nz_ind_gT_);
		}

		arg3_[0] = &xk_[0]; arg3_[1] = &uk_[0];
		for (int k = 0; k <= N_; k++){
			if (k != bookkeeper_->final_ind_){
				get_x(x, k); get_u(x, k);
				
				// fixed constraints
				if (blocks_->g_fixed_J_has_nz_){
					arg3_[2] = &p_g_fixed_[0];
					blocks_->g_fixed_J_(arg3_, res1_);
					addContribution(values,
									blocks_->g_fixed_J_.sparsity_out(0).nnz(),
									bookkeeper_->jac_nz_ind_g_fixed_[k]);
				}
			}
		}

		// extra constraints
		arg3_[0] = &xk_[0]; arg3_[1] = &uk_[0];
		for (const auto& pair : bookkeeper_->ind_g_e_){
			get_x(x, pair.first.k); get_u(x, pair.first.k);
			arg3_[2] = &p_g_extra_[pair.first.k][pair.first.j][pair.first.i][0];
			blocks_->g_extra_J_[pair.first.j](arg3_, res1_);
			addContribution(values, blocks_->g_extra_J_[pair.first.j].
												sparsity_out(0).nnz(),
							bookkeeper_->jac_nz_ind_g_extra_
								[pair.first.k][pair.first.j][pair.first.i]);
		}

		// dynamical constraints
		int k_init = 0;
		if (free_time_){
			arg5_[0] = &xx_(0,0); arg5_[1] = &uu_(0,0); arg5_[2] = &dts_[0];
			arg5_[3] = &x[0];
		} else {
			arg4_[0] = &xx_(0,0); arg4_[1] = &uu_(0,0); arg4_[2] = &dts_[0];
		}
		while (k_init != bookkeeper_->final_ind_){
			// cout<<"k = "<<k_init<<endl;
			if (free_time_){
				arg5_[4] = &p_g_disc_[bookkeeper_->
					ind_disc_block_[k_init].value()][0];
			} else {
				arg4_[3] = &p_g_disc_[bookkeeper_->
					ind_disc_block_[k_init].value()][0];
			}
		
			// get variable values
			getCollVars(k_init, bookkeeper_->nb_g_d_[k_init].value(), x);
			// get time-steps
			for (int i = 0; i < bookkeeper_->nb_g_d_[k_init].value()-1; i++){
				dts_[i] = bookkeeper_->dt(k_vals_[i]);
			}

			if (free_time_){
				blocks_->g_disc_J_[bookkeeper_->
					ind_disc_block_[k_init].value()](arg5_, res1_);
			} else {
				blocks_->g_disc_J_[bookkeeper_->
					ind_disc_block_[k_init].value()](arg4_, res1_);
			}

			addContribution(values, 
							blocks_->g_disc_J_
								[bookkeeper_->ind_disc_block_[k_init].value()]
								.sparsity_out(0).nnz(),
							bookkeeper_->jac_nz_ind_g_disc_[k_init]);
			k_init = k_vals_[bookkeeper_->nb_g_d_[k_init].value()-1];
		}
	}

	} catch (const std::exception &exc){std::cerr<<exc.what()<<endl; 
             std::exit(1);}
	if(print_){std::cout<<"exit "<<__func__<<std::endl;}
	return true;
};


bool NLPInterface::eval_h(
            Index         n,
            const Number* x,
            bool          new_x,
            Number        obj_factor,
            Index         m,
            const Number* lambda,
            bool          new_lambda,
            Index         nele_hess,
            Index*        iRow,
            Index*        jCol,
            Number*       values
){
	if(print_){std::cout<<"entering "<<__func__<<std::endl;}
	try{
	if (values == NULL){
		for (int i = 0; i < bookkeeper_->hess_nnz_; i++){
			iRow[i] = bookkeeper_->hess_row_ind_[i].value();
			jCol[i] = bookkeeper_->hess_col_ind_[i].value();
		}
	} else {
		if (!clearingStructuralZeros_){
			for (int i = 0; i < nele_hess; i++){values[i] = 0;}
		}
		if (free_time_){
			// initial cost
			arg3_[0] = &xk_[0]; arg3_[1] = &x[0]; arg3_[2] = &p_Phi_0_[0];
			get_x(x, 0);

			if (blocks_->Phi_0_H_has_nz_){
				blocks_->Phi_0_H_(arg3_, res1_);
				addContribution(values, blocks_->
									Phi_0_H_.sparsity_out(0).nnz(),
								bookkeeper_->hess_nz_ind_Phi_0_,
								blocks_->filter_Phi_0_H_, obj_factor);
			}

			// initial constraint
			if (blocks_->Phi_f_H_has_nz_){
				arg4_[0] = &xk_[0]; arg4_[1] = &x[0]; arg4_[2] = &lambda[0];
				arg4_[3] = &p_gT_[0];

				blocks_->g0_H_(arg4_, res1_);
				addContribution(values, blocks_->g0_H_.sparsity_out(0).nnz(),
								bookkeeper_->hess_nz_ind_g0_,
								blocks_->filter_g0_H_);
			}

			// final cost
			if (blocks_->Phi_f_H_has_nz_){
				arg3_[2] = &p_Phi_f_[0];
				get_x(x, bookkeeper_->final_ind_);

				blocks_->Phi_f_H_(arg3_, res1_);
				addContribution(values, blocks_->
									Phi_f_H_.sparsity_out(0).nnz(),
								bookkeeper_->hess_nz_ind_Phi_f_, 
								blocks_->filter_Phi_f_H_, obj_factor);
			}
		} else {
			// initial cost
			if (blocks_->Phi_0_H_has_nz_){
				arg2_[0] = &xk_[0]; arg2_[1] = &p_Phi_0_[0];
				get_x(x, 0);

				blocks_->Phi_0_H_(arg2_, res1_);
				addContribution(values, blocks_->
									Phi_0_H_.sparsity_out(0).nnz(),
								bookkeeper_->hess_nz_ind_Phi_0_, 
								blocks_->filter_Phi_0_H_, obj_factor);
			}

			// initial constraint
			if (blocks_->Phi_f_H_has_nz_){
				arg3_[0] = &xk_[0]; arg3_[1] = &lambda[0];
				arg3_[2] = &p_gT_[0];

				blocks_->g0_H_(arg3_, res1_);
				addContribution(values, blocks_->g0_H_.sparsity_out(0).nnz(),
								bookkeeper_->hess_nz_ind_g0_,
								blocks_->filter_g0_H_);
			}

			// final cost
			if (blocks_->Phi_f_H_has_nz_){
				arg2_[1] = &p_Phi_0_[0];
				get_x(x, bookkeeper_->final_ind_);

				blocks_->Phi_f_H_(arg2_, res1_);
				addContribution(values, blocks_->
									Phi_f_H_.sparsity_out(0).nnz(),
								bookkeeper_->hess_nz_ind_Phi_f_, 
								blocks_->filter_Phi_f_H_, obj_factor);
			}
		}

		// final constraint
		if (blocks_->gT_H_has_nz_){
			arg3_[0] = &xk_[0]; arg3_[1] = &lambda[blocks_->nb_g0_];
			arg3_[2] = &p_gT_[0];
			blocks_->gT_H_(arg3_, res1_);
			addContribution(values, blocks_->gT_H_.sparsity_out(0).nnz(),
							bookkeeper_->
								hess_nz_ind_gT_, blocks_->filter_gT_H_);
		}

		double dt;
		for (int k = 0; k <= N_; k++){
			if (k != bookkeeper_->final_ind_){
				get_x(x, k);
				get_u(x, k);
				arg4_[0] = &xk_[0]; arg4_[1] = &uk_[0]; 
				if (free_time_){
					arg5_[0] = &xk_[0]; arg5_[1] = &uk_[0];
					arg5_[2] = &dt; arg5_[3] = &x[0];
					arg5_[4] = &p_phi_[0]; 
				}

				// stage cost
				if (blocks_->phi_H_has_nz_){
					dt = bookkeeper_->dt(k);
					if (free_time_){
						blocks_->phi_H_(arg5_, res1_);
						addContribution(values,
										blocks_->phi_H_.sparsity_out(0).nnz(),
										bookkeeper_->hess_nz_ind_phi_[k], 
										blocks_->filter_phi_H_, obj_factor);
					} else {
						arg4_[2] = &dt; arg4_[3] = &p_phi_[0];
						blocks_->phi_H_(arg4_, res1_);
						addContribution(values,
										blocks_->phi_H_.sparsity_out(0).nnz(),
										bookkeeper_->hess_nz_ind_phi_[k], 
										blocks_->filter_phi_H_, obj_factor);
					}
				}

				// fixed constraint
				if (blocks_->g_fixed_H_has_nz_){
					arg4_[2] = &lambda[bookkeeper_->ind_g_f_[k].value()];
					arg4_[3] = &p_g_fixed_[0];
					blocks_->g_fixed_H_(arg4_, res1_);
					addContribution(values, 
									blocks_->g_fixed_H_.sparsity_out(0).nnz(),
									bookkeeper_->hess_nz_ind_g_fixed_[k],
									blocks_->filter_g_fixed_H_);
				}
			}
		}

		// extra constraint
		arg4_[0] = &xk_[0]; arg4_[1] = &uk_[0]; 
		for (const auto& pair : bookkeeper_->ind_g_e_){
			get_x(x, pair.first.k); get_u(x, pair.first.k);
			arg4_[2] = &lambda[pair.second];
			arg4_[3] = &p_g_extra_[pair.first.k][pair.first.j][pair.first.i][0];
			blocks_->g_extra_H_[pair.first.j](arg4_, res1_);
			addContribution(values, blocks_->g_extra_H_[pair.first.j].
												sparsity_out(0).nnz(),
							bookkeeper_->hess_nz_ind_g_extra_
								[pair.first.k][pair.first.j][pair.first.i],
							blocks_->filter_g_extra_H_[pair.first.j]);
		}

		// dynamical constraints
		int k_init = 0;

		if (free_time_){
			arg6_[0] = &xx_(0,0); arg6_[1] = &uu_(0,0); arg6_[2] = &dts_[0];
			arg6_[3] = &x[0];
		} else {
			arg5_[0] = &xx_(0,0); arg5_[1] = &uu_(0,0); arg5_[2] = &dts_[0];
		}
		while (k_init != bookkeeper_->final_ind_){
			if (free_time_){
				arg6_[5] = &p_g_disc_[bookkeeper_->
					ind_disc_block_[k_init].value()][0];
			} else {
				arg5_[4] = &p_g_disc_[bookkeeper_->
					ind_disc_block_[k_init].value()][0];
			}

			// get variable values
			getCollVars(k_init, bookkeeper_->nb_g_d_[k_init].value(), x);

			// get time-steps
			for (int i = 0; i < bookkeeper_->nb_g_d_[k_init].value()-1; i++){
				dts_[i] = bookkeeper_->dt(k_vals_[i]);
			}

			for (int i = 0; i < function_out_.size(); i++){
				function_out_[i] = -10;
			}

			if (free_time_){
				arg6_[4] = &lambda[bookkeeper_->ind_g_d_[k_init].value()];
				blocks_->g_disc_H_[bookkeeper_->
					ind_disc_block_[k_init].value()](arg6_, res1_);
			} else {
				arg5_[3] = &lambda[bookkeeper_->ind_g_d_[k_init].value()];
				blocks_->g_disc_H_[bookkeeper_->
					ind_disc_block_[k_init].value()](arg5_, res1_);
			}

			addContribution(values, 
							blocks_->g_disc_H_[bookkeeper_->
								ind_disc_block_[k_init].value()]
								.sparsity_out(0).nnz(),
							bookkeeper_->hess_nz_ind_g_disc_[k_init],
							blocks_->filter_g_disc_H_[bookkeeper_->
								ind_disc_block_[k_init].value()]);
			k_init = k_vals_[bookkeeper_->nb_g_d_[k_init].value()-1];
		}
	}

	} catch (const std::exception &exc){std::cerr<<exc.what()<<endl; 
             std::exit(1);}
	if(print_){std::cout<<"exit "<<__func__<<std::endl;}
	return true;
};

void NLPInterface::finalize_solution(
            SolverReturn               status,
            Index                      n,
            const Number*              x,
            const Number*              z_L,
            const Number*              z_U,
            Index                      m,
            const Number*              g,
            const Number*              lambda,
            Number                     obj_value,
            const IpoptData*           ip_data,
            IpoptCalculatedQuantities* ip_cq
){
	if(print_){std::cout<<"entering "<<__func__<<std::endl;}
	try{
	// store solution
	sol_obj_ = obj_value;
	for (int i = 0; i < n; i++){ sol_x_[i] = x[i];}
	for (int i = 0; i < m; i++){ sol_lambda_[i] = lambda[i];}

	} catch (const std::exception &exc){std::cerr<<exc.what()<<endl; 
             std::exit(1);}
	if(print_){std::cout<<"exit "<<__func__<<std::endl;}
};

void NLPInterface::updateN(int N){
	N_ = N;
	// any initial guess provided will no longer be valid
	initial_guess_available_ = false;
	sol_x_ = std::vector<double>(nx_*(N+1)+nu_*N+free_time_);
};

void NLPInterface::updateOrderedColsAndRows(){
};

void NLPInterface::updateParameters(std::vector<double>& p_Phi_0_,
									std::vector<double>& p_Phi_f_,
									std::vector<double>& p_phi_,
									std::vector<double>& p_g0_,
									std::vector<double>& p_gT_,
									std::vector<double>& p_g_fixed_){
	assert (blocks_->nb_p_Phi_0_ == 0 || 
			p_Phi_0_.size() == blocks_->nb_p_Phi_0_);
	assert (blocks_->nb_p_Phi_f_ == 0 || 
			p_Phi_f_.size() == blocks_->nb_p_Phi_f_);
	assert (blocks_->nb_p_phi_ == 0 || p_phi_.size() == blocks_->nb_p_phi_);
	assert (blocks_->nb_p_g0_ == 0 || p_g0_.size() == blocks_->nb_p_g0_);
	assert (blocks_->nb_p_gT_ == 0 || p_gT_.size() == blocks_->nb_p_gT_);
	assert (blocks_->nb_p_g_fixed_ == 0 || 
			p_g_fixed_.size() == blocks_->nb_p_g_fixed_);

	p_Phi_0_ = p_Phi_0_;
	p_Phi_f_ = p_Phi_f_;
	p_phi_ = p_phi_;
	p_g0_ = p_g0_;
	p_gT_ = p_gT_;
	p_g_fixed_ = p_g_fixed_;

	good_p_Phi_0_ = true;
	good_p_Phi_f_ = true;
	good_p_phi_ = true;
	good_p_g0_ = true;
	good_p_gT_ = true;
	good_p_g_fixed_ = true;
};

void NLPInterface::updateParameters(std::map<std::string, std::vector<double>>&
																parameters){
	bool invalid = false;
	for (const auto& pair : parameters) {
		if (pair.first == "p_Phi_0") {
			if (blocks_->nb_p_Phi_0_ == 0 || 
					pair.second.size() == blocks_->nb_p_Phi_0_){
				p_Phi_0_ = pair.second;
				good_p_Phi_0_ = true;
			} else {invalid = true; break;}
		} else if (pair.first == "p_Phi_f") {
			if (blocks_->nb_p_Phi_f_ == 0 || 
				pair.second.size() == blocks_->nb_p_Phi_f_){
			p_Phi_f_ = pair.second;
			good_p_Phi_f_ = true;
			} else {invalid = true; break;}
		} else if (pair.first == "p_phi") {
			if (blocks_->nb_p_phi_ == 0 || 
				   pair.second.size() == blocks_->nb_p_phi_){
			p_phi_ = pair.second;
			good_p_phi_ = true;
			} else {invalid = true; break;}
		} else if (pair.first == "p_g0") {
			if (blocks_->nb_p_g0_ == 0 || 
				   pair.second.size() == blocks_->nb_p_g0_){
			p_g0_ = pair.second;
			good_p_g0_ = true;
			} else {invalid = true; break;}
		} else if (pair.first == "p_gT") {
			if (blocks_->nb_p_gT_ == 0 || 
				   pair.second.size() == blocks_->nb_p_gT_){
			p_gT_ = pair.second;
			good_p_gT_ = true;
			} else {invalid = true; break;}
		} else if (pair.first == "p_g_fixed") {
			if (blocks_->nb_p_g_fixed_ == 0 || 
				pair.second.size() == blocks_->nb_p_g_fixed_){
			p_g_fixed_ = pair.second;
			good_p_g_fixed_ = true;
			} else {invalid = true; break;}
		}
	}
	if (invalid){
		throw illegalParameterProvided();
	}
};

void NLPInterface::updateDiscParameters(std::map<int, 
										std::vector<double>> parameters){
	for (const auto& pair : parameters){
		if (pair.first > blocks_->max_nb_steps_ || pair.first < 0 ||
				!blocks_->nb_steps_mapping_[pair.first].has_value()){
			throw unexistingParameterProvided();
		}
		if (blocks_->nb_p_g_disc_[pair.first] == 0 || 
				pair.second.size() == blocks_->nb_p_g_disc_[pair.first]){
			p_g_disc_[pair.first] = pair.second;
			good_p_g_disc_[pair.first] = true;
		}
	}
};

bool NLPInterface::checkParameters(){
	bool check = (good_p_Phi_0_ || blocks_->nb_p_Phi_0_ == 0) && 
				 (good_p_Phi_f_ || blocks_->nb_p_Phi_f_ == 0) && 
				 (good_p_phi_ || blocks_->nb_p_phi_ == 0) && 
				 (good_p_g0_ || blocks_->nb_p_g0_ == 0) && 
				 (good_p_gT_ || blocks_->nb_p_gT_ == 0) && 
				 (good_p_g_fixed_ || blocks_->nb_p_g_fixed_ == 0);
	
	for (int i = 0; i < blocks_->max_nb_steps_+1; i++){
		if (blocks_->nb_steps_mapping_[i].has_value()){
			check = check && 
				(blocks_->nb_p_g_disc_[blocks_->
					nb_steps_mapping_[i].value()] == 0 || good_p_g_disc_[i]);
		}
	}

	return check;
};

void NLPInterface::getCollVars(int k_init, int nk, const double* x, 
							   int k_break, int k_attach){
	get_x(x, k_init);
	get_u(x, k_init);
	xx_.col(0) = Eigen::Map<Eigen::VectorXd>(xk_.data(), xk_.size());
	uu_.col(0) = Eigen::Map<Eigen::VectorXd>(uk_.data(), uk_.size());
	k_vals_[0] = k_init;
	for (int j = 1; j < nk; j++){
		k_init = bookkeeper_->next(k_init, k_break, k_attach);
		k_vals_[j] = k_init;
		
		get_x(x, k_init);
		xx_.col(j) = Eigen::Map<Eigen::VectorXd>(xk_.data(), xk_.size());
		if (j < nk-1){
			get_u(x, k_init);
			uu_.col(j) = Eigen::Map<Eigen::VectorXd>(uk_.data(), uk_.size());
		}
	}
};

void NLPInterface::getCollVars(int k_init, int nk, const double* x){
	get_x(x, k_init);
	get_u(x, k_init);
	xx_.col(0) = Eigen::Map<Eigen::VectorXd>(xk_.data(), xk_.size());
	uu_.col(0) = Eigen::Map<Eigen::VectorXd>(uk_.data(), uk_.size());
	k_vals_[0] = k_init;
	for (int j = 1; j < nk; j++){
		k_init = bookkeeper_->next(k_init);
		k_vals_[j] = k_init;
		
		get_x(x, k_init);
		xx_.col(j) = Eigen::Map<Eigen::VectorXd>(xk_.data(), xk_.size());
		if (j < nk-1){
			get_u(x, k_init);
			uu_.col(j) = Eigen::Map<Eigen::VectorXd>(uk_.data(), uk_.size());
		}
	}
};

void NLPInterface::addContribution(double* values, int nnz,
                                   std::vector<std::optional<int>>& nz_ind,
								   std::vector<bool>& filter,
								   int multiplication_factor){
	int nz_ptr = 0;
	for (int i = 0; i < nnz; i++){
		if (filter[i]){ // ignore upper triangular values
			if (clearingStructuralZeros_ && 
					std::isnan(values[nz_ind[nz_ptr].value()])){
				values[nz_ind[nz_ptr].value()] = 0.0;
			}		
			values[nz_ind[nz_ptr].value()] += 
				multiplication_factor*function_out_[i];
			nz_ptr++;
		}
	}
};

void NLPInterface::clearStructuralZeros(){
	clearingStructuralZeros_ = true;

	clearVariableTraces();
	clearConstraintTraces();

	clearingStructuralZeros_ = false;
};

void NLPInterface::clearConstraintTraces(){
	/////////////////////////
	// make function calls //
	/////////////////////////
	int n = nx_*(N_+1) + nu_*N_ + free_time_;
	int m = bookkeeper_->next_free_ind_; 
	eval_g(n, &x_zeros_[0], true, m, &clear_constraints_[0]);
	eval_jac_g(n, &x_zeros_[0], true, m, bookkeeper_->jac_nnz_, nullptr, 
			   nullptr, &clear_jacobian_[0]);
	eval_h(n, &x_zeros_[0], true, 1.0, m, &lambda_zeros_[0], true, 
		   bookkeeper_->hess_nnz_, nullptr, nullptr, &clear_hessian_[0]);

	/////////////////////////////////////////////////////////////
	// compute mappings and shift vectors using these mappings //
	/////////////////////////////////////////////////////////////
	int nb_nz_cleared_g = initMapping(clear_constraints_, g_mapping_,
						  			  bookkeeper_->next_free_ind_);
	if (nb_nz_cleared_g > 0){
		shiftVectorElements(bookkeeper_->lb_container_, 
							bookkeeper_->next_free_ind_, g_mapping_);
		shiftVectorElements(bookkeeper_->ub_container_, 
							bookkeeper_->next_free_ind_, g_mapping_);
		bookkeeper_->next_free_ind_ -= nb_nz_cleared_g;
	}
	if (print_){
		std::cout<<"clearStructuralZeros(): cleared "<<nb_nz_cleared_g<<
			" constraint(s)"<<endl;
	}

	int nb_nz_cleared_jac = initMapping(clear_jacobian_, jac_mapping_,
						         		bookkeeper_->jac_nnz_);
	if (nb_nz_cleared_jac > 0){
		shiftVectorElements(bookkeeper_->jac_row_ind_, bookkeeper_->jac_nnz_,
							jac_mapping_);
		shiftVectorElements(bookkeeper_->jac_col_ind_, bookkeeper_->jac_nnz_,
							jac_mapping_);
		bookkeeper_->jac_nnz_ -= nb_nz_cleared_jac;
	}
	if (print_){
		std::cout<<"clearStructuralZeros(): cleared "<<nb_nz_cleared_jac<<
			" nonzero(s) in the jacobian"<<endl;
	}

	int nb_nz_cleared_hess = initMapping(clear_hessian_, hess_mapping_,
										 bookkeeper_->hess_nnz_);
	if (nb_nz_cleared_hess > 0){
		shiftVectorElements(bookkeeper_->hess_row_ind_, bookkeeper_->hess_nnz_,
							hess_mapping_);
		shiftVectorElements(bookkeeper_->hess_col_ind_, bookkeeper_->hess_nnz_,
							hess_mapping_);
		bookkeeper_->hess_nnz_ -= nb_nz_cleared_hess;
	}
	if (print_){
		std::cout<<"clearStructuralZeros(): cleared "<<nb_nz_cleared_hess<<
			" nonzero(s) in the hessian"<<endl;
	}
	
	///////////////////////////////////////////////
	// update pointer vectors using the mappings //
	///////////////////////////////////////////////
	if (nb_nz_cleared_g > 0){
		updateVectorElements(bookkeeper_->ind_g_f_, g_mapping_, N_+1);
		updateVectorElements(bookkeeper_->ind_g_d_, g_mapping_, N_+1);
		updateVectorElements(bookkeeper_->jac_row_ind_, g_mapping_,
							 bookkeeper_->jac_nnz_);
	}
	if (nb_nz_cleared_jac > 0){
		updateVectorElements(bookkeeper_->jac_nz_ind_g0_, jac_mapping_);
		updateVectorElements(bookkeeper_->jac_nz_ind_gT_, jac_mapping_);
	}
	if (nb_nz_cleared_hess > 0){
		updateVectorElements(bookkeeper_->hess_nz_ind_Phi_0_, hess_mapping_);
		updateVectorElements(bookkeeper_->hess_nz_ind_Phi_f_, hess_mapping_);
		updateVectorElements(bookkeeper_->hess_nz_ind_g0_, hess_mapping_);
		updateVectorElements(bookkeeper_->hess_nz_ind_gT_, hess_mapping_);
	}
	for (int k = 0; k <= N_; k++){
		if (k != bookkeeper_->final_ind_){
			if (nb_nz_cleared_jac > 0){
				updateVectorElements(bookkeeper_->jac_nz_ind_g_fixed_[k],
									jac_mapping_);
				if (k == bookkeeper_->ind_x_d_[k]){
					updateVectorElements(bookkeeper_->jac_nz_ind_g_disc_[k],
										 jac_mapping_);
				}
			}
			if (nb_nz_cleared_hess > 0){
				updateVectorElements(bookkeeper_->hess_nz_ind_phi_[k],
									hess_mapping_);
				updateVectorElements(bookkeeper_->hess_nz_ind_g_fixed_[k], 
									hess_mapping_);
				updateVectorElements(bookkeeper_->hess_nz_ind_g_disc_[k], 
									hess_mapping_);
			}
		}
	}
	for (auto& pair : bookkeeper_->ind_g_e_){
		if (nb_nz_cleared_g > 0){
			assert (g_mapping_[pair.second].has_value());
			pair.second = g_mapping_[pair.second].value();
		}
		if (nb_nz_cleared_jac > 0){
			updateVectorElements(bookkeeper_->
									jac_nz_ind_g_extra_[pair.first.k]
										[pair.first.j][pair.first.i],
									jac_mapping_);
		}
		if (nb_nz_cleared_hess > 0){
			updateVectorElements(bookkeeper_->
									hess_nz_ind_g_extra_[pair.first.k]
										[pair.first.j][pair.first.i],
									hess_mapping_);
		}
	}

	resetClearingVectors(nb_nz_cleared_g, nb_nz_cleared_jac,
						 nb_nz_cleared_hess);
}

void NLPInterface::clearVariableTraces(){
	if (!bookkeeper_->remove_traces_old_vars_){ return;}
	///////////////////
	// get k-mapping //
	///////////////////
	for (auto& e : k_mapping_){e.reset();}
	for (int i = 0; i < k_in_use_.size(); i++){
		k_in_use_[i] = std::numeric_limits<double>::quiet_NaN();
	}
	int k = 0;
	int max_k = bookkeeper_->final_ind_;
	while (k != bookkeeper_->final_ind_){
		max_k = std::max(max_k, k);
		k_in_use_[k] = 0.0;
		k = bookkeeper_->next_ind_[k].value();
	}
	k_in_use_[bookkeeper_->final_ind_] = 0.0;
	initMapping(k_in_use_, k_mapping_, Nmax_);
	int nb_k_cleared = getNbKCleared(k_mapping_);

	/////////////////////
	// get col mapping //
	/////////////////////
	for (int i = 0; i < col_in_use_.size(); i++){
		col_in_use_[i] = std::numeric_limits<double>::quiet_NaN();
	}
	k = 0;
	if (free_time_){
		col_in_use_[0] = 0.0;
	}
	while (k != bookkeeper_->final_ind_){
		std::fill(col_in_use_.begin() + (nx_+nu_)*k - 
						(k > bookkeeper_->final_ind_)*nu_ + free_time_,
				  col_in_use_.begin() + (nx_+nu_)*k - 
						(k > bookkeeper_->final_ind_)*nu_ + free_time_ + 
						nx_ + nu_,
				  0.0);
		k = bookkeeper_->next_ind_[k].value();
	}
	k = bookkeeper_->final_ind_;
	std::fill(col_in_use_.begin() + (nx_+nu_)*k + free_time_,
		      col_in_use_.begin() + (nx_+nu_)*k + free_time_ + nx_, 0.0);
	int nb_cols_cleared = initMapping(col_in_use_, col_mapping_, 
			   					      nx_*(max_k+1)+nu_*max_k+free_time_);

	///////////////////////////////////
	// shift vectors using k_mapping //
	///////////////////////////////////
	bool final_ind_changed = 
		bookkeeper_->final_ind_ != bookkeeper_->final_ind_old_;
	bookkeeper_->final_ind_ = k_mapping_[bookkeeper_->final_ind_].value();
	bookkeeper_->final_ind_old_ = bookkeeper_->final_ind_;
	shiftVectorElements(bookkeeper_->ind_g_f_, max_k+1, k_mapping_);
	shiftVectorElements(bookkeeper_->ind_g_d_, max_k+1, k_mapping_);
	shiftVectorElements(bookkeeper_->nb_g_d_, max_k+1, k_mapping_);
	shiftVectorElements(bookkeeper_->ind_disc_block_, max_k+1, k_mapping_);
	shiftVectorElements(bookkeeper_->ind_x_d_, max_k+1, k_mapping_);
	shiftVectorElements(bookkeeper_->next_ind_, max_k+1, k_mapping_);
	shiftVectorElements(bookkeeper_->time_from_ind_, max_k+1, k_mapping_);
	shiftVectorElements(bookkeeper_->jac_nz_ind_g_fixed_, max_k+1, k_mapping_);
	shiftVectorElements(bookkeeper_->jac_nz_ind_g_disc_, max_k+1, k_mapping_);
	shiftVectorElements(bookkeeper_->jac_nz_ind_g_extra_, max_k+1, k_mapping_);
	shiftVectorElements(bookkeeper_->hess_nz_ind_phi_, max_k+1, k_mapping_);
	shiftVectorElements(bookkeeper_->hess_nz_ind_g_fixed_, max_k+1, k_mapping_);
	shiftVectorElements(bookkeeper_->hess_nz_ind_g_disc_, max_k+1, k_mapping_);
	shiftVectorElements(bookkeeper_->hess_nz_ind_g_extra_, max_k+1, k_mapping_);
	shiftVectorElements(p_g_extra_, max_k+1, k_mapping_);

	if (final_ind_changed){
		bookkeeper_->ind_g_f_[bookkeeper_->final_ind_].reset();
		bookkeeper_->ind_g_d_[bookkeeper_->final_ind_].reset();
		bookkeeper_->nb_g_d_[bookkeeper_->final_ind_].reset();
		bookkeeper_->ind_disc_block_[bookkeeper_->final_ind_].reset();
		bookkeeper_->ind_x_d_[bookkeeper_->final_ind_].reset();
		bookkeeper_->next_ind_[bookkeeper_->final_ind_].reset();
		for (auto e : bookkeeper_->
				jac_nz_ind_g_fixed_[bookkeeper_->final_ind_]){
			e.reset();}
		for (auto e : bookkeeper_->
				jac_nz_ind_g_disc_[bookkeeper_->final_ind_]){
			e.reset();}
		for (auto e : bookkeeper_->
				hess_nz_ind_phi_[bookkeeper_->final_ind_]){
			e.reset();}
		for (auto e : bookkeeper_->
				hess_nz_ind_g_fixed_[bookkeeper_->final_ind_]){
			e.reset();}
		for (auto e : bookkeeper_->
				hess_nz_ind_g_disc_[bookkeeper_->final_ind_]){
			e.reset();}
	}
	
	//////////////////////////////
	// update extra constraints //
	//////////////////////////////
	auto it = bookkeeper_->ind_g_e_.begin();
	extraIndex updated_key;
	std::unordered_map<extraIndex, int> updated_map;
    while (it != bookkeeper_->ind_g_e_.end()) {
		updated_key = it->first;
        if (!k_mapping_[it->first.k].has_value() || 
				it->first.k == bookkeeper_->final_ind_) {
			bookkeeper_->nb_g_extra_applied_[it->first.j] -= 1;
            it = bookkeeper_->ind_g_e_.erase(it);
        } else if (it->first.k != k_mapping_[it->first.k].value()){
			updated_key.k = k_mapping_[it->first.k].value();
			if (updated_key.k != bookkeeper_->final_ind_){
				updated_map.insert(std::make_pair(updated_key, it->second));
			} else {
				bookkeeper_->nb_g_extra_applied_[it->first.j] -= 1;
			}
			it = bookkeeper_->ind_g_e_.erase(it);
        } else {
			++it;
		}
    }
	bookkeeper_->ind_g_e_.insert(updated_map.begin(), updated_map.end());

	///////////////////////////////////////////////
	// construct jac_mapping_ and shift elements //
	///////////////////////////////////////////////
	int nb_nz_cleared_jac = 0;
	for (int i = 0; i < bookkeeper_->jac_nnz_; i++){
		if (std::isnan(col_in_use_[bookkeeper_->jac_col_ind_[i].value()])){
			nb_nz_cleared_jac++;
			jac_mapping_[i].reset();
		} else {
			jac_mapping_[i] = i - nb_nz_cleared_jac;
		}
	}
	if (nb_nz_cleared_jac > 0){
		shiftVectorElements(bookkeeper_->jac_col_ind_, bookkeeper_->jac_nnz_, 
							jac_mapping_);
		shiftVectorElements(bookkeeper_->jac_row_ind_, bookkeeper_->jac_nnz_,
							jac_mapping_);
		bookkeeper_->jac_nnz_ -= nb_nz_cleared_jac;
	}

	////////////////////////////////////////////////
	// construct hess_mapping_ and shift elements //
	////////////////////////////////////////////////
	int nb_nz_cleared_hess = 0;
	for (int i = 0; i < bookkeeper_->hess_nnz_; i++){
		if (std::isnan(col_in_use_[bookkeeper_->hess_col_ind_[i].value()]) ||
				std::isnan(col_in_use_[bookkeeper_->hess_row_ind_[i].value()])){
			nb_nz_cleared_hess++;
			hess_mapping_[i].reset();
		} else {
			hess_mapping_[i] = i - nb_nz_cleared_hess;
		}
	}
	if (nb_nz_cleared_hess > 0){
		shiftVectorElements(bookkeeper_->hess_col_ind_, bookkeeper_->hess_nnz_, 
							hess_mapping_);
		shiftVectorElements(bookkeeper_->hess_row_ind_, bookkeeper_->hess_nnz_,
							hess_mapping_);
		bookkeeper_->hess_nnz_ -= nb_nz_cleared_hess;
	}
	
	/////////////////////
	// update contents //
	/////////////////////
	max_k -= nb_k_cleared;
	updateVectorElements(bookkeeper_->jac_col_ind_, col_mapping_,
						 bookkeeper_->jac_nnz_);
	updateVectorElements(bookkeeper_->hess_col_ind_, col_mapping_,
						 bookkeeper_->hess_nnz_);
	updateVectorElements(bookkeeper_->hess_row_ind_, col_mapping_,
						 bookkeeper_->hess_nnz_);
	updateVectorElements(bookkeeper_->ind_x_d_, k_mapping_, max_k+1);
	updateVectorElements(bookkeeper_->next_ind_, k_mapping_, max_k+1);
	if (nb_nz_cleared_jac > 0){
		updateVectorElements(bookkeeper_->jac_nz_ind_g0_, jac_mapping_);
		updateVectorElements(bookkeeper_->jac_nz_ind_gT_, jac_mapping_);
	}
	if (nb_nz_cleared_hess > 0){
		updateVectorElements(bookkeeper_->hess_nz_ind_Phi_0_, hess_mapping_);
		updateVectorElements(bookkeeper_->hess_nz_ind_Phi_f_, hess_mapping_);
		updateVectorElements(bookkeeper_->hess_nz_ind_g0_, hess_mapping_);
		updateVectorElements(bookkeeper_->hess_nz_ind_gT_, hess_mapping_);
	}

	for (int k = 0; k <= N_; k++){
		if (k != bookkeeper_->final_ind_){
			if (nb_nz_cleared_jac > 0){
				updateVectorElements(bookkeeper_->jac_nz_ind_g_fixed_[k],
									jac_mapping_);
				if (k == bookkeeper_->ind_x_d_[k]){
					updateVectorElements(bookkeeper_->jac_nz_ind_g_disc_[k],
										 jac_mapping_);
				}
			}
			if (nb_nz_cleared_hess > 0){
				updateVectorElements(bookkeeper_->hess_nz_ind_phi_[k],
									hess_mapping_);
				updateVectorElements(bookkeeper_->hess_nz_ind_g_fixed_[k], 
									hess_mapping_);
				updateVectorElements(bookkeeper_->hess_nz_ind_g_disc_[k], 
									hess_mapping_);
			}
		}
	}

	for (auto& pair : bookkeeper_->ind_g_e_){
		if (nb_nz_cleared_jac > 0){
			updateVectorElements(bookkeeper_->
									jac_nz_ind_g_extra_[pair.first.k]
										[pair.first.j][pair.first.i],
									jac_mapping_);
		}
		if (nb_nz_cleared_hess > 0){
			updateVectorElements(bookkeeper_->
									hess_nz_ind_g_extra_[pair.first.k]
										[pair.first.j][pair.first.i],
									hess_mapping_);
		}
	}

	bookkeeper_->remove_traces_old_vars_ = false;
};

void NLPInterface::resetClearingVectors(int nb_nz_cleared_g, 
										int nb_nz_cleared_jac, 
										int nb_nz_cleared_hess){
	for (int i = 0; i < bookkeeper_->next_free_ind_+nb_nz_cleared_g; i++){
		clear_constraints_[i] = std::numeric_limits<double>::quiet_NaN();
	}
	for (int i = 0; i < bookkeeper_->jac_nnz_+nb_nz_cleared_jac; i++){
		clear_jacobian_[i] = std::numeric_limits<double>::quiet_NaN();
	}
	for (int i = 0; i < bookkeeper_->hess_nnz_+nb_nz_cleared_hess; i++){
		clear_hessian_[i] = std::numeric_limits<double>::quiet_NaN();
	}

}

void NLPInterface::setInitialGuess(std::vector<double>& new_init_guess){
	assert (new_init_guess.size() == nx_*(N_+1) + nu_*N_ + free_time_);
	std::copy(new_init_guess.begin(), new_init_guess.end(),
			  init_guess_.begin());
	initial_guess_available_ = true;	
};

void NLPInterface::get_x(const double* x, int k){
	assert (k >= 0 && k <= N_);
	int offset = (k > bookkeeper_->final_ind_)*nu_;
	for (int i = 0; i < nx_; i++){
		xk_[i] = x[(nx_ + nu_)*k - offset + free_time_ + i];
	}
}

void NLPInterface::get_u(const double* x, int k){
	assert (k >= 0 && k <= N_);
	assert (k != bookkeeper_->final_ind_);
	int offset = (k > bookkeeper_->final_ind_)*nu_;
	for (int i = 0; i < nu_; i++){
		uk_[i] = x[(nx_ + nu_)*k + nx_ - offset + free_time_ + i];
	}
}

int NLPInterface::initMapping(std::vector<double>& v_out,
							  std::vector<std::optional<int>>& mapping,
							  int nnz){
	int nb_nz_cleared = 0;
	for (int i = 0; i < nnz; i++){
		if (std::isnan(v_out[i])){
			nb_nz_cleared++;
			mapping[i].reset();
		} else {
			mapping[i] = i - nb_nz_cleared;
		}
	}
	return nb_nz_cleared;
};

int NLPInterface::getNbKCleared(std::vector<std::optional<int>>& k_mapping_){
	int hole_counter = 0;
	int temp_hole_counter = 0;
	for (int k = 0; k < k_mapping_.size(); k++){
		if (k_mapping_[k].has_value()){
			hole_counter += temp_hole_counter;
			temp_hole_counter = 0;
		} else if (k != bookkeeper_->final_ind_){
			temp_hole_counter++;
		}
	}
	return hole_counter;
};

template <typename T>
void NLPInterface::shiftVectorElements(std::vector<std::optional<T>>& v,
									int nnz_old,
						 			std::vector<std::optional<int>>& mapping){
	for (int i = 0; i < nnz_old; i++){
		if (mapping[i].has_value()){
			v[mapping[i].value()] = v[i];
		}
	}
};

template <typename T>
void NLPInterface::shiftVectorElements(std::vector<T>& v,
									int nnz_old,
						 			std::vector<std::optional<int>>& mapping){
	for (int i = 0; i < nnz_old; i++){
		if (mapping[i].has_value()){
			v[mapping[i].value()] = v[i];
		}
	}
};

template <typename T>
void NLPInterface::updateVectorElements(std::vector<std::optional<T>>& v,
									std::vector<std::optional<int>>& mapping,
									int max_ind){
	for (int i = 0; i < max_ind; i++){
		if (v[i].has_value()){
			// If this assertion fails, it means vector v has a value pointing
			// somewhere that does not exist anymore
			assert (mapping[v[i].value()].has_value());
			v[i] = mapping[v[i].value()].value();
		}
	}
};
template <typename T>
void NLPInterface::updateVectorElements(std::vector<std::optional<T>>& v,
									std::vector<std::optional<int>>& mapping){
	for (int i = 0; i < v.size(); i++){
		if (v[i].has_value()){
			assert (mapping[v[i].value()].has_value());
			v[i] = mapping[v[i].value()].value();
		}
	}
};


void NLPInterface::copyFunctionOutToTarget(double* target, int target_ind,
										   std::vector<int>& mapping,
										   bool includes_time){
	int offset = 0;
	if (includes_time && mapping[0] == 0){
		target[0] += function_out_[0];
		offset = 1;
	}
	for (int i = offset; i < mapping.size(); i++){
		target[target_ind + mapping[i]] += function_out_[i];
	}
}
