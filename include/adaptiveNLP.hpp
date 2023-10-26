#ifndef __ADAPTIVE_NLP__
#define __ADAPTIVE_NLP__

#include "../include/buildingBlocks.hpp"
#include "../include/bookkeeper.hpp"
#include "../include/NLPInterface.hpp"
#include "IpIpoptApplication.hpp"
#include <casadi/casadi.hpp>
#include <vector>
#include <optional>
#include <numeric>
#include <iostream>

#ifndef NDEBUG
#include "../include/interfaceTester.hpp"
#endif

using namespace casadi;

// Class to provide functionanlity to change NLP problem structures at low
// cost. 
class AdaptiveNLP{
    public:
        AdaptiveNLP()
        :Nmax_(0), nx_(0), nu_(0), free_time_(false), blocks_(nullptr) {};

        // constructor
        // blocks:                  building blocks to be used to assemble the 
        //                          NLP
        // T:                       horizon lenght. For free-time problems, it
        //                          should be equal to 1.0
        // Nmax_:                   maximum number of time-steps to be 
        //                          accounted for in the NLP
        // max_nb_extra_instances_: maximum number of times the same extra 
        //                          constraint can be applied to the same
        //                          time-step
        AdaptiveNLP(BuildingBlocks& blocks, double T, int Nmax, 
                    int max_nb_extra_instances_);

        ~AdaptiveNLP(){};

        // Function to initialize some timesteps.
        // The length of nks is the number of intervals to be added to the
        // problem. nks[i] specifies how many steps are involved in interval i
        void initTimeSteps(std::vector<int> nks, std::vector<double> tt = {});

        // Function adds time-steps to the current sequence in the same 
        // horizon (points are squezed in by updating the time stamps of all 
        // time-steps).
        // The length of nks is the number of intervals to be added to the
        // problem. nks[i] specifies how many steps are involved in interval i
        void addTimeSteps(std::vector<int> nks, std::vector<double> tt = {});

        // Function extends the horizon by appending time-steps to the current
        // sequence.
        // The length of nks is the number of intervals to be added to the
        // problem. nks[i] specifies how many steps are involved in interval i
        void appendTimeSteps(std::vector<int> nks,
                             std::vector<double> tt = {});

        // Function to reduce the horizon length by a specified number of i
        // discretization intervals
        void reduceHorizon(int nb_intervals);

        // Function to change the discretization functions used. The interval
        // starting at time-step k is replaced by a number of intervals equal
        // to the length of nks. The number of time-steps involved in each of 
        // these new intervals is determined by nks[i]
        void changeIntervalDiscretization(int k, std::vector<int> nks,
                                          std::vector<double> tt = {});

        // Function adds the indicated extra constraints on the given indices.
        // If no constraint indices are provided, zero indices are used.
        // If no parameters are provided, the default parameter (0) is used
        // If a vector of length 1 is provided for the constraint_ind and the
        // parameters, this constraint index and parameter vector is used for
        // all indices to which to apply the constraint.
        void addExtraConstraint(std::vector<int> inds,
                                std::vector<int> constraint_ind,
                                std::vector<int> instance_ind,
                                std::vector<std::vector<double>> parameters);
        void addExtraConstraint(std::vector<int> inds,
                                std::vector<int> constraint_ind,
                                std::vector<std::vector<double>> parameters){
            std::vector<int> instance_ind(constraint_ind.size(), -1);
            return addExtraConstraint(inds, constraint_ind, instance_ind, 
                                      parameters);
        };
        void addExtraConstraint(std::vector<int> inds){
            std::vector<int> constraint_ind(inds.size(), 0);
            std::vector<std::vector<double>> parameters(inds.size(), {0});
            return addExtraConstraint(inds, constraint_ind, parameters);
        };
        void addExtraConstraint(std::vector<int> inds,
                                std::vector<std::vector<double>> parameters){
            std::vector<int> constraint_ind(parameters.size(), 0);
            return addExtraConstraint(inds, constraint_ind, parameters);
        };
        void addExtraConstraint(std::vector<int> inds,
                                std::vector<int> constraint_ind){
            std::vector<std::vector<double>> parameters(constraint_ind.size(),
                                                        {0});
            return addExtraConstraint(inds, constraint_ind, parameters);
        };

        // Function to remove extra constraints. If no instance indices are
        // provided, the first instance is removed. If no constraint index is
        // provided either, the first extra constraint is removed.
        // If an instance index (or constraint index) of -1 is provided, the 
        // first intance (or constraint) is removed.
        void removeExtraConstraint(std::vector<int> inds){
            std::vector<int> constraint_ind = std::vector<int>(inds.size(),-1);
            std::vector<int> instance_ind = std::vector<int>(inds.size(),-1);
            return removeExtraConstraint(inds, constraint_ind, instance_ind);
        };
        void removeExtraConstraint(std::vector<int> inds,
                                  std::vector<int> constraint_ind){
            std::vector<int> instance_ind = std::vector<int>(inds.size(),-1);
            return removeExtraConstraint(inds, constraint_ind, instance_ind);
        };
        void removeExtraConstraint(std::vector<int> inds,
                                   std::vector<int> constraint_ind,
                                   std::vector<int> instance_ind);

        // Function to remove extra constraints. If no instance index is
        // provided, the first instance is removed. If no constraint index is
        // provided either, the first extra constraint is removed.
        // If an instance index (or constraint index) of -1 is provided, the 
        // first intance (or constraint) is removed.
        void removeExtraConstraint(std::vector<int> inds, int constraint_ind){
            return removeExtraConstraint(inds, constraint_ind, -1);
        };
        void removeExtraConstraint(std::vector<int> inds, int constraint_ind,
                                   int instance_ind);

        // Function to solve the NLP. The relevant parameters have to be
        // provided.
        // If a correct output variable is provided, the computation time can
        // be written in it.
        // An initial guess can also be specified.
        int solveNlp(std::map<std::string, std::vector<double>> parameters){
            double computation_time;
            return solveNlp(parameters, computation_time);
        };
        int solveNlp(std::map<std::string, std::vector<double>> parameters,
                                     std::vector<double>& init_guess){
            double computation_time;
            return solveNlp(parameters, computation_time, init_guess);
        };
        int solveNlp(std::map<std::string, std::vector<double>> parameters,
                     double& computation_time, 
                     std::vector<double>& init_guess){
            interface_->setInitialGuess(init_guess);
            return solveNlp(parameters, computation_time);
        };
        int solveNlp(std::map<std::string, std::vector<double>> parameters, 
                     double& computation_time);

        std::vector<double> getSolution();

        #ifndef NDEBUG
        void testInterface(){
            std::vector<double> x = std::vector<double>(nx_*(N_+1) + nu_*N_ + 
                                                            free_time_, 1.0);
            std::cout<<
                "[AdaptiveNLP::testInterface()] getting initial guess"<<endl;
            interface_->initial_guess_available_ = false;
            interfaceTester_.get_starting_point(x);
            testInterface(x);
        };
        void testInterface(std::vector<double>& x);

        std::vector<double> computeFunctionEvaluationTimes(){
            return interfaceTester_.computeFunctionEvaluationTimes();
        };

        #endif

        // Function to clear the structural zeros in the interface
        void clearStructuralZeros(){ interface_->clearStructuralZeros();};

        // basic getters
        int getNx(){ return nx_;};
        int getNu(){ return nu_;};
        int getN(){ return N_;};
        int getNmax(){ return Nmax_;};
        bool getFreeTime(){ return free_time_;};
        double getT(){ return T_;};
        int getNumberOfConstraints(){return bookkeeper_.next_free_ind_;};
        int getNumberOfBoundaryConstraints(){
            return blocks_->nb_g0_ + blocks_->nb_g0_;
        };
        int getNumberOfFixedConstraints(){
            return blocks_->nb_g_fixed_*N_;
        };
        int getNumberOfDiscretizationConstraints(){
            return blocks_->nx_*N_;
        };
        std::vector<int> getNumberOfExtraConstraints(){
            return bookkeeper_.nb_g_extra_applied_;
        };
        std::vector<int> getJacRows(){
            std::vector<int> res(bookkeeper_.jac_nnz_);
            for (int i = 0; i < bookkeeper_.jac_nnz_; i++){
                res[i] = bookkeeper_.jac_row_ind_[i].value();
            }
            return res;
        }
        std::vector<int> getJacCols(){
            std::vector<int> res(bookkeeper_.jac_nnz_);
            for (int i = 0; i < bookkeeper_.jac_nnz_; i++){
                res[i] = bookkeeper_.jac_col_ind_[i].value();
            }
            return res;
        }
        std::vector<int> getHessRows(){
            std::vector<int> res(bookkeeper_.hess_nnz_);
            for (int i = 0; i < bookkeeper_.hess_nnz_; i++){
                res[i] = bookkeeper_.hess_row_ind_[i].value();
            }
            return res;
        }
        std::vector<int> getHessCols(){
            std::vector<int> res(bookkeeper_.hess_nnz_);
            for (int i = 0; i < bookkeeper_.hess_nnz_; i++){
                res[i] = bookkeeper_.hess_col_ind_[i].value();
            }
            return res;
        }
        int getFinalInd(){ return bookkeeper_.final_ind_;};
        std::vector<std::optional<int>> getNextInd(){
            std::vector<std::optional<int>> 
                next_ind(bookkeeper_.next_ind_.size());
            std::copy(bookkeeper_.next_ind_.begin(),
                      bookkeeper_.next_ind_.end(), next_ind.begin());
            return next_ind;
        }
        int getNext(int k){
            if (bookkeeper_.next_ind_[k].has_value()){
                return bookkeeper_.next_ind_[k].value();
            } else {return -1;}
        }
        std::vector<std::optional<double>> getTimeFromIndex(){
            std::vector<std::optional<double>> 
                time_from_ind(bookkeeper_.time_from_ind_.size());
            std::copy(bookkeeper_.time_from_ind_.begin(),
                      bookkeeper_.time_from_ind_.end(),
                      time_from_ind.begin());
            return time_from_ind;
        }
        double getTimeFromIndex(int k){
            assert(bookkeeper_.time_from_ind_[k].has_value());
            return bookkeeper_.time_from_ind_[k].value();
        }
        int getNbSteps(int k){
            if (bookkeeper_.nb_g_d_[k].has_value()){
                return bookkeeper_.nb_g_d_[k].value();
            } else {return -1;}
        }
        std::vector<std::optional<int>> getNbSteps(){
            return bookkeeper_.nb_g_d_;
        }
        double getIntervalLength(int k){
            assert (bookkeeper_.nb_g_d_[k].has_value());
            double res = 0.0;
            for (int i = 0; i < bookkeeper_.nb_g_d_[k].value()-1; i++){
                res += bookkeeper_.dt(k);
                if (i < bookkeeper_.nb_g_d_[k].value()-2){
                    k = bookkeeper_.next(k);
                }
            }
            return res;
        }
        int getIterCount(){ return interface_->getIterCount();};
        double getNanoSecondCounter(){ return interface_->getNanoSecondCounter();}
        
        void setIpoptPrintLevel(int level){print_level_ = level;}

    private:
        // Function to check input to functions that change the problem 
        // structure (such as adding variables)
        // illegalChangeToNLP exceptions will be thrown in cases:
        //      - the problem is not initialized yet
        //      - the provided list of number of time-steps is empty
        //      - there is an interval with invalid number of time-steps
        //      - there is not enough space allocated to include the new steps
        // This function also returns the nb of time-steps to add
        int checkInputToNLPUpdate(std::vector<int>& nks);

        // Function adds the stage cost and fixed constraints to the number of 
        // stages provided starting at index k_start
        // It is assumed that the bookkeeping (next_ind) is already udpated
        void addStageContributions(int k_start, int nb_of_stages);

        // Function adds the stage cost and fixed constraints to the number of 
        // stages provided starting at index k_start
        // It is assumed that the bookkeeping (next_ind and time_from_index)
        // is already updated 
        void addDiscretizationContributions(int k_start, 
                                            std::vector<int>& nks);

        // Function to return the indices involved in the discretization
        // constraint starting at k_init and involving nk time-steps.
        // The result is written in class attribute 'k_vals_'
        void getIndices(int k_init, int nk);

        // Function to return the row and column indices of the nonzero values
        // of the give sparsity. 
        // If 'only_lower_triangular' = true, only the lower triangular 
        // elements are returned.
        // Optional arguments row_offset and col_offset can be used to locate
        // the given sparisty in a larger matrix, such that a correct decision
        // can be made to find lower triangular elements
        // The resulting row and column indices (indices local to this 
        // sparsity) are written in class attributes 'nz_rows_' and 'nz_cols_'
        int getNonzeroIndices(const Sparsity& sp, bool only_lower_triangular){
            return getNonzeroIndices(sp, only_lower_triangular, 0, 0);
        };
        int getNonzeroIndices(const Sparsity& sp, bool only_lower_triangular,
                              int row_offset, int col_offset);

        // Function to add an elementary sparsity contribution to the sparsity
        // defined by row_ind and col_ind.
        // An elementary contribution is one that is compact (meaning this
        // sparsity can be added as is without having to split it)
        int addElementaryContribution(Sparsity& sp, int row_offset,
                                    int col_offset, bool hessian, 
                                    std::vector<std::optional<int>>& row_ind,
                                    std::vector<std::optional<int>>& col_ind);

        // Function to add the given sparsity contribution to the jacobian
        // sparsity.
        // sp:                sparsity to add
        // k_indices:         time-steps involved in this contribution
        // includes_time:     does this contribution have a time contribution?
        // row_offset:        the first row-index of this contribution
        // bookkeeper_vector: vector in which indices will be written such
        //                    that the NLPInterface object can easily put
        //                    numerical values in the correct place in the
        //                    values list.
        void addCompleteJacobianContribution(const Sparsity& sp, 
                                             std::vector<int> k_indices, 
                                             bool includes_time,
                                             int row_offset,
                                             std::vector<std::optional<int>>&
                                                            bookkeeper_vector);

        // Function to add the given sparsity contribution to the hessian
        // sparsity.
        // sp:                sparsity to add
        // k_indices:         time-steps involved in this contribution
        // includes_time:     does this contribution have a time contribution?
        // row_offset:        the first row-index of this contribution
        // bookkeeper_vector: vector in which indices will be written such
        //                    that the NLPInterface object can easily put
        //                    numerical values in the correct place in the
        //                    values list.
        //                    The i-th element in the bookkeeper_vector
        //                    contains the sparsity index of the i-th
        //                    lower triangular nonzero element of the sparsity
        //                    where the nonzeros are traversed in a column-
        //                    major manner. 
        void addCompleteHessianContribution(const Sparsity& sp,
                                            std::vector<int> k_indices,
                                            bool includes_time,
                                            std::vector<std::optional<int>>& 
                                                            bookkeeper_vector);

        // Function to compute the maximal number of constraints to account
        // for in the NLP
        int getMaxNbConstraints(int max_nb_extra_instances_){
            return blocks_->nb_g0_ + blocks_->nb_gT_ + 
                   (Nmax_)*(blocks_->nb_g_fixed_ + 
                            std::accumulate(blocks_->nb_g_extra_.begin(),
                                            blocks_->nb_g_extra_.end(), 0)*
                                            max_nb_extra_instances_ +
                            nx_);
        }

        // Function finds the first instance of any extra constraint at the 
        // given time-step. It returns whether it has found one. If it did,
        // the result is written in the given variables 'constraint_ind' and
        // 'instance_ind'
        bool getFirstConstraintAndInstanceInd(int k, int& constraint_ind,
                                              int& instance_ind);

        // Function finds the first instance of the given extra constraint at
        // the given time-step. It returns whether it has found one. If it did,
        // the result is written in the given variable 'instance_ind'
        bool getFirstInstanceInd(int k, int constraint_ind, int& instance_ind);

        // Function returns the first free instance slot for the specified
        // constraint at the given time-step k. It returns whether it has
        // found one. If it did, the result is written in the given variable
        // 'instance_ind'
        bool getFreeInstanceSlot(int k, int constraint_ind, int& instance_ind);

        // current number of time-steps in the problem
        int N_;

        // maximal number of time-steps to account for
        const int Nmax_;

        // number of states
        const int nx_;

        // number of inputs
        const int nu_;

        // length of the time horizon
        double T_;

        // boolean indicating whether the problem is free-time or not
        const bool free_time_;

        // pointer to the building blocks object (always points to the same
        // object)
        BuildingBlocks* const blocks_;

        // Bookkeeper object
        Bookkeeper bookkeeper_;

        // NLPInterface interface;
        SmartPtr<IpoptApplication> ipopt_app_;
        SmartPtr<NLPInterface> interface_;
        int print_level_ = 5;

        ///////////////////////////////
        // scratch space allocations //
        ///////////////////////////////
        std::vector<int> k_vals_;
        std::vector<int> nz_rows_;
        std::vector<int> nz_cols_;
        std::vector<int> nz_indices_;

        // for addCompleteHessianContribution
        // rows to be selected for the subsparsity
        std::vector<casadi_int> selected_rows_nxnu_;
        std::vector<casadi_int> selected_rows_nx_;
        std::vector<casadi_int>* row_ptr_;

        // columns to be selected for the subsparsity
        std::vector<casadi_int> selected_cols_nxnu_;
        std::vector<casadi_int> selected_cols_nx_;
        std::vector<casadi_int>* col_ptr_;
        
        //empty mapping for the subsparsity
        std::vector<casadi_int> mapping_ = std::vector<casadi_int>(0);

        #ifndef NDEBUG
        InterfaceTester interfaceTester_;
        #endif
};

#endif