#ifndef __BOOKKEEPER__
#define __BOOKKEEPER__
#include <vector>
#include <optional>
#include "buildingBlocks.hpp"
#include <unordered_map>
#include <tuple>

struct extraIndex{
    int k; int j; int i;

    bool operator==(const extraIndex& other) const{
        return k == other.k && j == other.j && i == other.i;
    }
};
namespace std {
    template<>
    struct hash<extraIndex> {
        size_t operator()(const extraIndex index) const {
            return hash<int>()(index.k) ^ hash<int>()(index.j) ^ 
                hash<int>()(index.i);
        }
    };
};

class NLPInterface;
class AdaptiveNLP;

// Class to perform collect bookkeeping information. This information is
// updated by the AdaptiveNLP class and used by the NLPInterface class.
class Bookkeeper{
    friend class NLPInterface;
    friend class AdaptiveNLP;

    public:
        // default constructor
        Bookkeeper():
        Nmax_(0), nx_(0), nu_(0), max_nb_constraints_(0), nb_extra_blocks_(0),
        nb_disc_blocks_(0), max_nb_extra_instances_(0), free_time_(false){};

        // constructor
        // blocks_:                 collection of buildling blocks
        // Nmax_  :                 limit on the total number of time steps in
        //                          the NLP to be accounted for
        // max_nb_constraints_:     total number of constraints in the NLP to 
        //                          be accounted for
        // max_nb_extra_instances_: total amount of time the same extra
        //                          constraints can be applied to the same 
        //                          time-step (e.g. with different parameters) 
        Bookkeeper(BuildingBlocks& blocks, int Nmax, int max_nb_constraints,
                   int max_nb_extra_instances);

        // returns the time between x_{k} and x_{k+1}. If T is provided, the
        // time-grid is multiplied by this value
        double dt(int k, double T = 1.0);

        // return the index of the next time-step chronologically
        int next(int k, int k_break, int k_attach);
        int next(int k);

        // return the index of the previous time-step chronologically
        int previous(int k);

    protected:
        // lower and upper bounds on constraints
        std::vector<double> lb_container_;
        std::vector<double> ub_container_;

        // current number of extra constraints
        std::vector<int> nb_g_extra_applied_;

        // jacobian sparsities
        std::vector<std::optional<int>> jac_row_ind_;
        std::vector<std::optional<int>> jac_col_ind_;
        int jac_nnz_;

        // hessian sparsities
        std::vector<std::optional<int>> hess_row_ind_;
        std::vector<std::optional<int>> hess_col_ind_;
        int hess_nnz_;

        // indices in the jacobian sparsity of nonzero elements of different
        // contributions to the jacobian
        std::vector<std::optional<int>> jac_nz_ind_g0_;
        std::vector<std::optional<int>> jac_nz_ind_gT_;
        std::vector<std::vector<std::optional<int>>> jac_nz_ind_g_fixed_;
        std::vector<std::vector<std::optional<int>>> jac_nz_ind_g_disc_;
        std::vector<std::vector<std::vector<std::vector<std::optional<int>>>>>
                                                        jac_nz_ind_g_extra_;
                                    
        // indices in the hessian sparsity of nonzero elements of different
        // contributions to the hessian
        std::vector<std::optional<int>> hess_nz_ind_Phi_0_;
        std::vector<std::optional<int>> hess_nz_ind_Phi_f_;
        std::vector<std::vector<std::optional<int>>> hess_nz_ind_phi_;
        std::vector<std::optional<int>> hess_nz_ind_g0_;
        std::vector<std::optional<int>> hess_nz_ind_gT_;
        std::vector<std::vector<std::optional<int>>> hess_nz_ind_g_fixed_;
        std::vector<std::vector<std::optional<int>>> hess_nz_ind_g_disc_;
        std::vector<std::vector<std::vector<std::vector<std::optional<int>>>>>
                                                        hess_nz_ind_g_extra_;

        // indices in the constraint vector of fixed constraints
        std::vector<std::optional<int>> ind_g_f_;

        // indices in the constraint vector of extra constraints
        std::unordered_map<extraIndex, int> ind_g_e_;

        // indices in the constraint vector of discretization constraints
        std::vector<std::optional<int>> ind_g_d_;
        // number of time-steps involved in the discretization constraint
        std::vector<std::optional<int>> nb_g_d_;
        // index of the discretization constraint corresponding to the correct
        // number of time-steps chosen
        std::vector<std::optional<int>> ind_disc_block_;
        // index of the first time-step involved in the discretization
        // constraint
        std::vector<std::optional<int>> ind_x_d_;

        // next free index in the constraint vector. is also used to derive
        // the current number of constraints in the NLP
        int next_free_ind_;

        // index of the (chronologically) next index
        std::vector<std::optional<int>> next_ind_;

        // index of the final state
        int final_ind_;
        int final_ind_old_; // to be used when reducing the horizon

        // timestamp of state xk
        std::vector<std::optional<double>> time_from_ind_;

        // maximum number of times the same extra constraint can be applied
        // to the same (xk, uk)-pair
        const int max_nb_extra_instances_;

        // maximum number of constraints accounted for in total
        const int max_nb_constraints_;

        bool remove_traces_old_vars_ = false;
        bool remove_traces_old_constraints_ = false;

    private:
        // total number of time steps to be accounted for
        const int Nmax_;
        
        // number of states
        const int nx_;
        
        // number of inputs
        const int nu_;
        
        // number of extra constraints that exist
        const int nb_extra_blocks_;

        // number of discretization blocks
        const int nb_disc_blocks_;

        // boolean indicating if this is a free-time problem or not
        const bool free_time_;
};

#endif