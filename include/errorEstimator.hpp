#ifndef __ERROR_ESTIMATOR__
#define __ERROR_ESTIMATOR__
#include <cassert>
#include <vector>
#include <casadi/casadi.hpp>

using namespace casadi;

// Object to estimate the relative error in the system dynamics when direct 
// collocation is used.
class errorEstimator{
    public:
        // Constructor
        errorEstimator();

        // Constructor
        // nx is the number of states, nu is the number of controls and fd is 
        // a casadi function to evaluate the dynamics that hase three input 
        // args (x, u, p) and has nx outputs.
        errorEstimator(int nx_, int nu_, Function& fd_);

        void updateVars(int nx_, int nu_, Function fd_);

        // Function that returns the relative error estimate on the system 
        // dynamics.
        // The vector tau indicates the collocation points mapped on the 
        // interval (0, 1].
        // xx is a vector where every element is a vector representing the 
        // states at a (non)collocated point. uu is defined in a similar 
        // way but for the controls.
        // The value tf provided here is the length (in time) of this single 
        // interval.
        // Also for free-time problems, tf provided is the true length in time.
        double getEstimate(double tf_, std::vector<double> tau, 
                           std::vector<std::vector<double>>& xx,
                           std::vector<std::vector<double>>& uu);

        // Function to evaluate the solution on the collocation points tau 
        // represented by xx and uu on the new collocation points given by 
        // tauHat. The solution is written in xxHat and uuHat.
        void evaluateNewGrid(std::vector<double>& tau, 
							 std::vector<double>& tauHat,
                			 std::vector<std::vector<double>>& xx, 
							 std::vector<std::vector<double>>& xxHat,
                			 std::vector<std::vector<double>>& uu, 
							 std::vector<std::vector<double>>& uuHat);

        // Function to construct an integration matrix. The result is stored 
		// in a class attribute.
        void constructIntegrationMatrix(std::vector<double>& tauHat);

    private:

        // Function to construct an interpolating polynomial for the states 
		// and the inputs.
        void constructInterpolant(std::vector<double>& tau, 
								  std::vector<std::vector<double>>& xx,
                                  std::vector<std::vector<double>>& uu);

        // Function to retrieve the integrated values that are used for error 
		// estimation. The resulting values yy[i][j] represent the values 
		// int(p_F(t), t = t0...tau_i) for state j where p_F represents the 
		// polynomial through the evaluations of the dynamics at the newly 
		// sampled points
        void getIntegratedValues(std::vector<double>& tauHat, 
								 std::vector<std::vector<double>>& xxHat,
                				 std::vector<std::vector<double>>& uuHat, 
								 std::vector<std::vector<double>>& yy);

        // number of states
        int nx;

        // number of inputs
        int nu;

        // vector of interpolating polynomials for the states
        std::vector<Polynomial> px;

        // vector of interpolating polynomials for the controls
        std::vector<Polynomial> pu;

        // Integration matrix where B_i,j = int(lagrange_j(t), t = t0...tau_i)
        std::vector<std::vector<double>> B;

        // Function to evaluate the system dynamics with inputs (x, u, p) and 
		// nx outputs
        Function fd;

        // Lenght of the time interval
        double tf;
};

#endif