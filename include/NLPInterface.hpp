#ifndef __NLPINTERFACE__
#define __NLPINTERFACE__

#include "buildingBlocks.hpp"
#include "bookkeeper.hpp"
#include <IpTNLP.hpp>
#include <cassert>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <map>

using namespace Ipopt;

class AdaptiveNLP;

// Class to provide the interface with the solver (IPOPT)
class NLPInterface : public TNLP{
    friend class AdaptiveNLP;

    public:
        // Default constructor
        NLPInterface():
        print_(false), nx_(0), nu_(0), free_time_(false),
        Nmax_(0), blocks_(NULL), bookkeeper_(NULL){};

        // Constructor
        // N:              current number of timesteps in the NLP
        // bookkeeper:     bookkeeper object
        // buildingBlocks: collection of building blocks to assemble the NLP
        // diagonalize:    whether or not to diagonalize the NLP 
        //                 (TODO: use this)
        NLPInterface(int N, Bookkeeper& bookkeeper,
                     BuildingBlocks& buildingBlocks, 
                     bool diagonalize);
        NLPInterface(int N, Bookkeeper& bookkeeper,
                     BuildingBlocks& buildingBlocks)
            :NLPInterface(N, bookkeeper, buildingBlocks, false)
            {};

        ~NLPInterface(){};

        // Function to provide basic info about the NLP to be solved
        virtual bool get_nlp_info(
            Index&          n,
            Index&          m,
            Index&          nnz_jac_g,
            Index&          nnz_h_lag,
            IndexStyleEnum& index_style
       );

        // Function to provide bounds of constraints and on variables
        virtual bool get_bounds_info(
            Index   n,
            Number* x_l,
            Number* x_u,
            Index   m,
            Number* g_l,
            Number* g_u
        );

        // Function to provide initial guess
        virtual bool get_starting_point(
            Index   n,
            bool    init_x,
            Number* x,
            bool    init_z,
            Number* z_L,
            Number* z_U,
            Index   m,
            bool    init_lambda,
            Number* lambda
        );

        // Function to evaluate the objective function
        virtual bool eval_f(
            Index         n,
            const Number* x,
            bool          new_x,
            Number&       obj_value
        );

        // Function to evaluate the objective function gradient
        virtual bool eval_grad_f(
            Index         n,
            const Number* x,
            bool          new_x,
            Number*       grad_f
        );

        // Function to evaluate the constraints
        virtual bool eval_g(
            Index         n,
            const Number* x,
            bool          new_x,
            Index         m,
            Number*       g
        );

        // Function to evaluate the constraint jacobian
        virtual bool eval_jac_g(
            Index         n,
            const Number* x,
            bool          new_x,
            Index         m,
            Index         nele_jac,
            Index*        iRow,
            Index*        jCol,
            Number*       values
        );

        // Function to evaluate the lagrangian hessian
        virtual bool eval_h(
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
        );

        // Function to finalize the solution
        virtual void finalize_solution(
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
        );

        virtual bool intermediate_callback(
            AlgorithmMode              mode,
            Index                      iter,
            Number                     obj_value,
            Number                     inf_pr,
            Number                     inf_du,
            Number                     mu,
            Number                     d_norm,
            Number                     regularization_size,
            Number                     alpha_du,
            Number                     alpha_pr,
            Index                      ls_trials,
            const IpoptData*           ip_data,
            IpoptCalculatedQuantities* ip_cq
        ){
            iter_count_ = iter;
            return true;
        };

        int getIterCount(){ return iter_count_;};
        double getNanoSecondCounter(){ return nanosecond_counter;}

    protected:
        // Function to change the number of time-steps in the NLP
        void updateN(int N);

        // TODO
        void updateOrderedColsAndRows();

        // Function to update the parameters of the different functions of the
        // NLP
        void updateParameters(std::vector<double>& p_Phi_0_,
                              std::vector<double>& p_Phi_f_,
                              std::vector<double>& p_phi_,
                              std::vector<double>& p_g0_,
                              std::vector<double>& p_gT_,
                              std::vector<double>& p_g_fixed_);

        // Function to update the parameters of the different functions of the
        // NLP. The dictonary keys contain the name of the parameter and the
        // values contain the new parameter value.
        // Possible names are 'p_Phi_0', 'p_Phi_f', 'p_phi', 'p_g0', 'p_gT',
        // and 'p_g_fixed'
        void updateParameters(std::map<std::string, std::vector<double>>&
                                parameters);

        // Function to update the parameters of the discretization functions.
        // The dictionary contains an integer indicating the number of 
        // timesteps involved in the constraint of which the parameters need
        // updating and the vector contains the parameter values.
        void updateDiscParameters(std::map<int,
                                  std::vector<double>> parameters);

        // Fucntion that checks whether all relevant parameters have been
        // initialized with correct sizes.
        bool checkParameters();

        // Function to return the state vectors and control vectors involved 
        // in the discretization constraint starting with state x_{k_init}
        // that contains the specified number of time-steps involved.
        // The k-indices are written to the attritbute 'k_vals_'
        void getCollVars(int k_init, int nk, const double* x,
                         int k_break, int k_attach);
        void getCollVars(int k_init, int nk, const double* x);

        // Function to add a contribution to the vector of nonzero elements
        // values: vector of nonzero values to be updated
        // nnz:    number of nonzero elements to be written
        // nz_ind: vector that contains the index in the vector 'values' for 
        //         every nonzero value to be written. This vector is created
        //         by AdaptiveNLP
        // filter: vector indicating for every nonzero element whether it
        //         should be added or not. It is used for hessian 
        //         contributions to only add lower-triangular values. For
        //         jacobian contributions, the 'filter_default_' attribute of 
        //         the blocks_ object can be used (which is true for all
        //         nonzero elements).
        void addContribution(double* values, int nnz,
                             std::vector<std::optional<int>>& nz_ind,
                             int multiplication_factor=1){
            addContribution(values, nnz, nz_ind, blocks_->filter_default_,
                            multiplication_factor);
        };
        void addContribution(double* values, int nnz,
                             std::vector<std::optional<int>>& nz_ind,
                             std::vector<bool>& filter,
                             int multiplication_factor=1);

        // Function to remove constraints that have been removed in the past
        void clearStructuralZeros();

        // Function to remove constraint traces in the bookkeeping
        void clearConstraintTraces();

        // Function to update all bookkeeping to account for removed variables 
        void clearVariableTraces();

        // Function to set an initial guess
        void setInitialGuess(std::vector<double>& newinit_guess);

    private:
        // Function to retrieve the vector xk from the complete solution
        // vector.
        // The result is written in class attribute xk_
        void get_x(const double* x, int k);

        // Function to retrieve the vector uk from the complete solution
        // vector
        // The result is written in class attribute uk_
        void get_u(const double* x, int k);

        void resetClearingVectors(int nb_nz_cleared_g=0, 
                                  int nb_nz_cleared_jac=0, 
                                  int nb_nz_cleared_hess=0);

        // Function to initialize the shifting mappings in
        // clearStructuralZeros(). For the i-th element in v_out, mapping[i]
        // contains the index in which this element should be placed to stack
        // all relevant values to the front/
        // The function returns the number of "holes" have been found in v_out
        int initMapping(std::vector<double>& v_out,
                        std::vector<std::optional<int>>& mapping,
                        int nnz);

        // Function to compute how many 'holes' there are in the k_mapping
        // because of removing variables. If variables at the end are removed,
        // these are not counted.
        // This function will give a different result compared to initMapping
        // because it is unknown how many variables have been removed so no
        // proper 'nnz' can be provided to initMapping, making the resulting
        // hole count unreliable.
        int getNbKCleared(std::vector<std::optional<int>>& k_mapping_);

        // Function to shift elements of a vector such that all elements that
        // have a value are stacked to the front of the value using the
        // provided mapping.
        template <typename T>
        void shiftVectorElements(std::vector<std::optional<T>>& v,
                                 std::vector<std::optional<int>>& mapping){
            return shiftVectorElements(v, v.size(), mapping);
        };
        template <typename T>
        void shiftVectorElements(std::vector<std::optional<T>>& v, int nnz_old,
                                 std::vector<std::optional<int>>& mapping);
        template <typename T>
        void shiftVectorElements(std::vector<T>& v, int nnz_old,
                                 std::vector<std::optional<int>>& mapping);

        // Function to update the content of the given vectors such that they
        // point to the updated indices using the given mapping. Instead of 
        // poinint to i, the content of v will point to mapping[i]
        template <typename T>
        void updateVectorElements(std::vector<std::optional<T>>& v,
                                  std::vector<std::optional<int>>& mapping,
                                  int max_ind);
        template <typename T>
        void updateVectorElements(std::vector<std::optional<T>>& v,
                                  std::vector<std::optional<int>>& mapping);

        // Function to copy the contents of class attribute 'function_out_' to
        // the target vector, starting at index target_ind. 
        // NOTE: The target_ind is the index without considering the fact that
        // the problem might be free-time!
        // For every nonzero of the function output vector, 'mapping' provides
        // the index in the output vector.
        // If 'includes_time'=true, the first element of 'function_out_' is
        // always written to target[0].
        void copyFunctionOutToTarget(double* target, int target_ind, 
                                     std::vector<int>& mapping, 
                                     bool includes_time=false);

        // variable to enable or disable printing information
        // const bool print_;
        bool print_;
    
        const int nx_;
        const int nu_;
        const bool free_time_;
        int N_;
        const int Nmax_;
        int iter_count_;

        bool initial_guess_available_ = false;
        std::vector<double> init_guess_;

        // pointer to collection of building blocks
        BuildingBlocks* const blocks_;

        // pointer to bookkeeper object
        Bookkeeper* const bookkeeper_;
        
        // boolean indicating whether to diagonalize the hessian. If set to
        // true, a permutation will be done on the variable vector of IPOPT
        // and on the column indices of the jacobian and on the column and
        // row indices of the hessian.
        bool diagonalize_;

        ///////////////////////////////
        // scratch space allocations //
        ///////////////////////////////
        std::vector<int> k_vals_;    // k_vals_ from getCollVars
        std::vector<double> dts_;    // list of time indexes
        std::vector<double> xk_;     // state vector
        std::vector<double> uk_;     // input vector
        Eigen::MatrixXd xx_;         // used in getCollVars
        Eigen::MatrixXd uu_;         // used in getCollVars
        // row indices of nonzero elements of a contribution
        std::vector<int> nz_rows_;
        // column indices of nonzero elements of a contribution
        std::vector<int> nz_cols_;

        // input argument pointers for casadi functions
        std::vector<const double*> arg1_ = std::vector<const double*>(1);
        std::vector<const double*> arg2_ = std::vector<const double*>(2);
        std::vector<const double*> arg3_ = std::vector<const double*>(3);
        std::vector<const double*> arg4_ = std::vector<const double*>(4);
        std::vector<const double*> arg5_ = std::vector<const double*>(5);
        std::vector<const double*> arg6_ = std::vector<const double*>(6);
        
        // allocated space for casadi functions output
        std::vector<double> function_out_;
        
        // output argument pointer for casadi functions
        const std::vector<double*> res1_;

        // scratch space for "clear structural zeros"
        std::vector<double> x_zeros_;
        std::vector<double> lambda_zeros_;
        std::vector<double> clear_constraints_;
        std::vector<double> clear_jacobian_;
        std::vector<double> clear_hessian_;
        std::vector<std::optional<int>> g_mapping_;
        std::vector<std::optional<int>> jac_mapping_;
        std::vector<std::optional<int>> hess_mapping_;
        bool clearingStructuralZeros_ = false;

        std::vector<double> k_in_use_;
        std::vector<std::optional<int>> k_mapping_;
        std::vector<double> col_in_use_;
        std::vector<std::optional<int>> col_mapping_;

        // function parameters
        std::vector<double> p_Phi_0_;
        bool good_p_Phi_0_ = false;
        std::vector<double> p_Phi_f_;
        bool good_p_Phi_f_ = false;
        std::vector<double> p_phi_;
        bool good_p_phi_ = false;
        std::vector<double> p_g0_;
        bool good_p_g0_ = false;
        std::vector<double> p_gT_;
        bool good_p_gT_ = false;
        std::vector<double> p_g_fixed_;
        bool good_p_g_fixed_ = false;
        std::vector<std::vector<double>> p_g_disc_;
        std::vector<bool> good_p_g_disc_;
        std::vector<std::vector<std::vector<std::vector<double>>>> p_g_extra_;

        // solution storage
        std::vector<double> sol_x_;
        std::vector<double> sol_lambda_;
        double sol_obj_;

        double nanosecond_counter = 0.0;
};

#endif