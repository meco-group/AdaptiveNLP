#ifndef __BUILDING_BLOCKS__
#define __BUILDING_BLOCKS__

#include <casadi/casadi.hpp>
#include <optional>
#include <cassert>

using namespace casadi;

class NLPInterface;
class Bookkeeper;
class AdaptiveNLP;
class AdaptiveCorridorHelper;
class ManyObstaclesHelper;
class VirtualObstacleHelper;
class MoonlanderHelper;

// Class to collect the building blocks defining a discretized OCP.
// For fixed-time problems, the discrete NLP is given by:
//
//      min   Phi_0(x0) + sum_{0}^{N-1} phi(xk, uk, dt) + Phi_f(xN)
//       s.t.      gDisc(xk, uk, dt, xk+1) == 0
//                 lb0 <= g0(x0) <= ub0
//                 lbT <= gT(xN) <= ubT
//                 lbf <= gf(xk,uk) <= ubf
//                 lbe <= ge(xk,uk) <= ube
//
//             (all functions also include a parameter input)
//
// For free-time problems, the disctrete NLP is given by:
//
//      min   Phi_0(x0, T) + sum_{0}^{N-1} phi(xk, uk, dt) + Phi_f(xN, T)
//       s.t.      gDisc(xk, uk, dt, xk+1) == 0
//                 lb0 <= g0(x0, T) <= ub0
//                 lbT <= gT(xN) <= ubT
//                 lbf <= gf(xk,uk) <= ubf
//                 lbe <= ge(xk,uk) <= ube
//
//             (all functions also include a parameter input)
class BuildingBlocks{
    friend class NLPInterface;
    friend class Bookkeeper;
    friend class AdaptiveNLP;
    friend class AdaptiveCorridorHelper;
    friend class ManyObstaclesHelper;
    friend class VirtualObstacleHelper;
    friend class MoonlanderHelper;

    public:
        // default Constructor
        BuildingBlocks();

        // Add problem sizes
        // nx is the number of states, nu is the number of controls, T is the
        // length of the time horizon. For free-time problems, T is typically
        // equal to one (scaled time)
        void add_sizes(int nx, int nu, bool free_time=false);

        // Add stage cost
        // The function signatures are given by:
        // free-time:
        // phi:   in:  [xk (nx x 1), uk (nu x 1), dt (1 x 1), p (nb_p x 1)]
        //        out: 1 x 1
        // phi_J: in:  [xk (nx x 1), uk (nu x 1), dt (1 x 1), p (nb_p x 1)]
        //        out: 1 x (1 + nx + nu)
        // phi_H: in:  [xk (nx x 1), uk (nu x 1), dt (1 x 1), p (nb_p x 1)]
        //        out: (1 + nx + nu) x (1 + nx + nu)
        // (dt is scaled by the total-time)
        //
        // fixed-time:
        // phi:   in:  [xk (nx x 1), uk (nu x 1), dt (1 x 1), p (nb_p x 1)]
        //        out: 1 x 1
        // phi_J: in:  [xk (nx x 1), uk (nu x 1), dt (1 x 1), p (nb_p x 1)]
        //        out: 1 x (nx + nu)
        // phi_H: in:  [xk (nx x 1), uk (nu x 1), dt (1 x 1), p (nb_p x 1)]
        //        out: (nx + nu) x (nx + nu)
        //
        // The jacobian and hessian ordering of columns is given by
        // T - xk - uk (or xk - uk)
        void add_phi(Function& phi, Function& phi_J, Function& phi_H,
                     int nb_p=0);

        // Add cost to initial state
        // The function signatures are given by:
        // free-time:
        // Phi_0:   in:  [xk (nx x 1), T (1 x 1), p (nb_p x 1)]
        //          out: 1 x 1
        // Phi_0_J: in:  [xk (nx x 1), T (1 x 1), p (nb_p x 1)]
        //          out: 1 x (1 + nx)
        // Phi_0_H: in:  [xk (nx x 1), T (1 x 1), p (nb_p x 1)]
        //          out: (1 + nx) x (1 + nx)
        //
        // fixed-time:
        // Phi_0:   in:  [xk (nx x 1), p (nb_p x 1)]
        //          out: 1 x 1
        // Phi_0_J: in:  [xk (nx x 1), p (nb_p x 1)]
        //          out: 1 x nx
        // Phi_0_H: in:  [xk (nx x 1), p (nb_p x 1)]
        //          out: nx x nx
        //
        // The jacobian and hessian ordering of columns is given by
        // T - xk (or xk)
        void add_Phi_0(Function& Phi_0, Function& Phi_0_J, Function& Phi_0_H,
                       int nb_p=0);

        // Add cost to final state
        // The function signatures are given by:
        // free-time:
        // Phi_f:   in:  [xk (nx x 1), T (1 x 1), p (nb_p x 1)]
        //          out: 1 x 1
        // Phi_f_J: in:  [xk (nx x 1), T (1 x 1), p (nb_p x 1)]
        //          out: 1 x (1 + nx)
        // Phi_f_H: in:  [xk (nx x 1), T (1 x 1), p (nb_p x 1)]
        //          out: (1 + nx) x (1 + nx)
        //
        // fixed-time:
        // Phi_f:   in:  [xk (nx x 1), p (nb_p x 1)]
        //          out: 1 x 1
        // Phi_f_J: in:  [xk (nx x 1), p (nb_p x 1)]
        //          out: 1 x nx
        // Phi_f_H: in:  [xk (nx x 1), p (nb_p x 1)]
        //          out: nx x nx
        //
        // The jacobian and hessian ordering of columns is given by
        // T - xk (or xk)
        void add_Phi_f(Function& Phi_f, Function& Phi_f_J, Function& Phi_f_H,
                       int nb_p=0);

        // Add initial conditions
        // The function signatures are given by:
        // free-time:
        // g0:   in:  [xk (nx x 1), T (1 x 1), p (nb_p x 1)]
        //       out: m x 1
        // g0_J: in:  [xk (nx x 1), T (1 x 1), p (nb_p x 1)]
        //       out: m x (1 + nx)
        // g0_H: in:  [xk (nx x 1), T (1 x 1), lambda (m x 1), p (nb_p x 1)]
        //       out: (1 + nx) x (1 + nx)
        //
        // fixed-time:
        // g0:   in:  [xk (nx x 1), p (nb_p x 1)]
        //       out: m x 1
        // g0_J: in:  [xk (nx x 1), p (nb_p x 1)]
        //       out: m x nx
        // g0_H: in:  [xk (nx x 1), lambda (m x 1), p (nb_p x 1)]
        //       out: nx x nx
        //
        // The jacobian and hessian ordering of columns is given by
        // T - xk (or xk)
        void add_g0(Function& g0, Function& g0_J, Function& g0_H,
                    std::vector<double>& lb, std::vector<double>& ub,
                    int nb_p=0);

        // Add final conditions
        // The function signatures are given by:
        //
        // gT:   in:  [xk (nx x 1), p (nb_p x 1)]
        //       out: m x 1
        // gT_J: in:  [xk (nx x 1), p (nb_p x 1)]
        //       out: m x nx
        // gT_H: in:  [xk (nx x 1), lambda (m x 1), p (nb_p x 1)]
        //       out: nx x nx
        //
        // The jacobian and hessian ordering of columns is given by xk
        void add_gT(Function& gT, Function& gT_J, Function& gT_H,
                    std::vector<double>& lb, std::vector<double>& ub,
                    int nb_p=0);

        // Add fixed constraints
        // The function signatures are given by:
        //
        // g_fixed:   in:  [xk (nx x 1), uk (nu x 1), p (nb_p x 1)]
        //            out: m x 1
        // g_fixed_J: in:  [xk (nx x 1), uk (nu x 1), p (nb_p x 1)]
        //            out: m x (nx + nu)
        // g_fixed_H: in:  [xk (nx x 1), uk (nu x 1), lambda (m x 1), 
        //                  p (nb_p x 1)]
        //            out: (nx + nu) x (nx + nu)
        //
        // The jacobian and hessian ordering of columns is given by
        // xk - uk
        void add_g_fixed(Function& g_fixed, Function& g_fixed_J,
                         Function& g_fixed_H, std::vector<double>& lb,
                         std::vector<double>& ub, int nb_p=0);

        // Add extra constraints
        // Lists of functions are provided. For every triplet of functions
        // g_extra[i], g_extra_J[i], g_extra_H[i], the function signatures are
        // given by:
        //
        // g_extra:   in:  [xk (nx x 1), uk (nu x 1), p (nb_p x 1)]
        //            out: m x 1
        // g_extra_J: in:  [xk (nx x 1), uk (nu x 1), p (nb_p x 1)]
        //            out: m x (nx + nu)
        // g_extra_H: in:  [xk (nx x 1), uk (nu x 1), lambda (m x 1), 
        //                  p (nb_p x 1)]
        //            out: (nx + nu) x (nx + nu)
        //
        // The jacobian and hessian ordering of columns is given by
        // xk - uk
        void add_g_extra(std::vector<Function> g_extra,
                         std::vector<Function> g_extra_J,
                         std::vector<Function> g_extra_H,
                         std::vector<std::vector<double>> lb,
                         std::vector<std::vector<double>> ub,
                         std::vector<int> nb_p){
            assert (nb_p.size() == g_extra.size());
            add_g_extra(g_extra, g_extra_J, g_extra_H, lb, ub);
            nb_p_g_extra_ = nb_p;
        }
        void add_g_extra(std::vector<Function> g_extra,
                         std::vector<Function> g_extra_J,
                         std::vector<Function> g_extra_H,
                         std::vector<std::vector<double>> lb,
                         std::vector<std::vector<double>> ub);
                        
        // Add discretization constraints
        // Lists of functions are provided. For every triplet of functions
        // g_disc[i], g_disc_J[i], g_disc_H[i], the function signatures are
        // given by:
        // The function signatures are given by:
        // free-time:
        // g_disc:   in:  [xx (nx x n+1), uu (nu x n), dt (n x 1), T (1 x 1),
        //                 p (nb_p x 1)]
        //           out: m x 1
        // g_disc_J: in:  [xx (nx x n+1), uu (nu x n), dt (n x 1), T (1 x 1),
        //                 p (nb_p x 1)]
        //           out: m x (1 + nx*(n+1) + nu*n)
        // g_disc_H: in:  [xx (nx x n+1), uu (nu x n), dt (n x 1), T (1 x 1),
        //                 lambda (m x 1), p (nb_p x 1)]
        //           out: (1 + nx*(n+1) + nu*n) x (1 + nx*(n+1) + nu*n)
        //
        // fixed-time:
        // g_disc:   in:  [xx (nx x n+1), uu (nu x n), dt (n x 1),
        //                 p (nb_p x 1)]
        //           out: m x 1
        // g_disc_J: in:  [xx (nx x n+1), uu (nu x n), dt (n x 1),
        //                 p (nb_p x 1)]
        //           out: m x (nx*(n+1) + nu*n)
        // g_disc_H: in:  [xx (nx x n+1), uu (nu x n), dt (n x 1),
        //                 lambda (m x 1), p (nb_p x 1)]
        //           out: (nx*(n+1) + nu*n) x (nx*(n+1) + nu*n)
        //
        // The jacobian and hessian ordering of columns is given by
        // T - x(k) - u(k) - x(k+1) - u(k+1) ... x(m+1) 
        // (or x(k) - u(k) - x(k+1) - u(k+1) ... x(m+1))
        void add_g_disc(std::vector<Function> g_disc, 
                        std::vector<Function> g_disc_J,
                        std::vector<Function> g_disc_H, 
                        std::vector<int> nb_p){
            assert (nb_p.size() == g_disc.size());
            add_g_disc(g_disc, g_disc_J, g_disc_H);
            nb_p_g_disc_ = nb_p;
        };
        void add_g_disc(std::vector<Function> g_disc, 
                        std::vector<Function> g_disc_J,
                        std::vector<Function> g_disc_H);

        // Add initial guesses
        // The function signatures are given by:
        // x_init: in:  [t (1 x 1)]
        //         out: nx x 1
        // u_init: int: [t (1 x 1)]
        //         out: nu x 1
        void add_inits(Function& x_init, Function& u_init);
        void add_inits(Function& x_init, Function& u_init, double t_init);

        // Function checks whether all building blocks have been initialised
        // and returns true it that is the case
        bool check_complete();

        // Function returns true if the given sizes match the values of this
        // collection of building blocks
        bool check_problem_sizes(int nx_check, int nu_check);

        // get nx
        int get_nx(){ return nx_;};

        // get nu
        int get_nu(){ return nu_;};

        // get free-time
        bool getFreeTime(){ return free_time_;};

    protected:
        // problem sizes
        int nx_;
        int nu_;

        // function returns the number of nonzero elements in the output of 
        // a jacobian. The input 'block_nb' specifies the jacobian.
        // 0: g0_J
        // 1: gT_J
        // 2: g_fixed_J
        int get_nb_jac_nnz(int block_nb);

        // function returns the number of nonzero elements in the output of 
        // the jacobian of an extra constraint. The input 'extra_nb'
        // specifies the index of the extra constraint.
        int get_nb_jac_nnz_extra(int extra_nb);

        // function returns the number of nonzero elements in the output of 
        // the jacobian of a discretization constraint. The input 'disc_nb'
        // specifies the index of the discretization constraint.
        int get_nb_jac_nnz_disc(int disc_nb);

        // function returns the number of nonzero elements in the output of 
        // a hessian. The input 'block_nb' specifies the hessian.
        // 0: Phi_0_H
        // 1: phi_H
        // 2: Phi_f_H
        // 3: g0_H
        // 4: gT_H
        // 5: g_fixed_H
        int get_nb_hess_nnz(int block_nb);

        // function returns the number of nonzero elements in the output of 
        // the hessian of an extra constraint. The input 'extra_nb'
        // specifies the index of the extra constraint.
        int get_nb_hess_nnz_extra(int extra_nb);

        // function returns the number of nonzero elements in the output of 
        // the hessian of a discretization constraint. The input 'disc_nb'
        // specifies the index of the discretization constraint.
        int get_nb_hess_nnz_disc(int disc_nb);

        bool phi_J_has_nz_;
        bool phi_H_has_nz_;

        bool Phi_0_J_has_nz_;
        bool Phi_0_H_has_nz_;

        bool Phi_f_J_has_nz_;
        bool Phi_f_H_has_nz_;

        // number of initial constraints
        int nb_g0_;
        bool g0_J_has_nz_;
        bool g0_H_has_nz_;

        // number of final constraints
        int nb_gT_;
        bool gT_J_has_nz_;
        bool gT_H_has_nz_;

        // number of fixed constraints
        int nb_g_fixed_;
        bool g_fixed_J_has_nz_;
        bool g_fixed_H_has_nz_;

        // list of numbers of extra constraints
        std::vector<int> nb_g_extra_;
        int nb_g_extra_total_ = 0;
        std::vector<bool> g_extra_J_has_nz_;
        bool g_extra_J_has_nz_total_;
        std::vector<bool> g_extra_H_has_nz_;
        bool g_extra_H_has_nz_total_;

        // number of extra constraints
        int nb_extra_blocks_;

        // number of discretization constraints
        std::vector<int> nb_g_disc_;
        std::vector<bool> g_disc_J_has_nz_;
        bool g_disc_J_has_nz_total_;
        std::vector<bool> g_disc_H_has_nz_;
        bool g_disc_H_has_nz_total_;
        int max_disc_J_nnz_;
        int max_disc_H_nnz_;

        // number of disc blocks
        int nb_disc_blocks_;

        // max number of time-steps involved in a discretization constraint
        int max_nb_steps_ = 0;

        // mapping that maps a number of steps to a discretization function
        // used as: g_disc_n = g_disc_[nb_steps_mapping_[n]]
        std::vector<std::optional<int>> nb_steps_mapping_;

        // max number of constraints over all constraint blocks
        int max_nb_g_ = 0;

        // numbers of parameters
        int nb_p_phi_;
        int nb_p_Phi_0_;
        int nb_p_Phi_f_;
        int nb_p_g0_;
        int nb_p_gT_;
        int nb_p_g_fixed_;
        std::vector<int> nb_p_g_disc_;
        std::vector<int> nb_p_g_extra_;

        // largest number of nonzeros in one function
        int max_nnz_ = 0;

        // mapppings from nonzeros to indices (for each nonzero, what index 
        // does it have in the output vector of the function?)
        std::vector<int> mapping_Phi_0_J_;
        std::vector<int> mapping_Phi_f_J_;
        std::vector<int> mapping_phi_J_;

        // flags to indicate if a nonzero element is lower triangular (and
        // hence should be added to the hessian) or not
        std::vector<bool> filter_Phi_0_H_;
        std::vector<bool> filter_Phi_f_H_;
        std::vector<bool> filter_phi_H_;
        std::vector<bool> filter_g0_H_;
        std::vector<bool> filter_gT_H_;
        std::vector<bool> filter_g_fixed_H_;
        std::vector<std::vector<bool>> filter_g_disc_H_;
        std::vector<std::vector<bool>> filter_g_extra_H_;

        // default filter with all true values
        std::vector<bool> filter_default_;
 
    private:
        // initialize the mapping vectors
        void initMappingVector(const Sparsity& sp, std::vector<int>& mapping);

        // initialize the filter vector
        void initFilterVector(const Sparsity& sp, std::vector<bool>& filter);

        // problem parameters
        bool free_time_;
        bool sizes_initialized_ = false;

        // stage cost
        Function phi_;
        Function phi_J_;
        Function phi_H_;
        bool stage_cost_initialized_ = false;

        // initial cost
        Function Phi_0_;
        Function Phi_0_J_;
        Function Phi_0_H_;
        bool initial_cost_initialized_ = false;

        // final cost
        Function Phi_f_;
        Function Phi_f_J_;
        Function Phi_f_H_;
        bool final_cost_initialized_ = false;

        // initial constraint
        Function g0_;
        Function g0_J_;
        Function g0_H_;
        std::vector<double> g0_lb_;
        std::vector<double> g0_ub_;
        bool initial_constraint_initialized_ = false;

        // final constraint
        Function gT_;
        Function gT_J_;
        Function gT_H_;
        std::vector<double> gT_lb_;
        std::vector<double> gT_ub_;
        bool final_constraint_initialized_ = false;

        // fixed constraint
        Function g_fixed_;
        Function g_fixed_J_;
        Function g_fixed_H_;
        std::vector<double> g_fixed_lb_;
        std::vector<double> g_fixed_ub_;
        bool fixed_constraint_initialized_ = false;

        // extra constraints
        std::vector<Function> g_extra_;
        std::vector<Function> g_extra_J_;
        std::vector<Function> g_extra_H_;
        std::vector<std::vector<double>> g_extra_lb_;
        std::vector<std::vector<double>> g_extra_ub_;
        bool extra_constraint_initialized_ = false;

        // disctretization constraints
        std::vector<Function> g_disc_;
        std::vector<Function> g_disc_J_;
        std::vector<Function> g_disc_H_;
        bool discretization_constraint_initialized_ = false;

        // initialization functions
        Function x_init_;
        Function u_init_;
        double t_init_;
        bool init_functions_initialized_ = false;
};

#endif