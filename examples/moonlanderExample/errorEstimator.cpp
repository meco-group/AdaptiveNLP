#include "errorEstimator.hpp"
#include <vector>
#include <casadi/casadi.hpp>
#include <cstdlib>

using namespace casadi;
using namespace std;

// Constructor
errorEstimator::errorEstimator(){
    nx = 1;
    nu = 1;
    tf = 1.0;
    fd = Function("fd", {MX::sym("x",1,1)}, {0.0});
    px = std::vector<Polynomial>(nx);
    pu = std::vector<Polynomial>(nu);
};

errorEstimator::errorEstimator(int nx_, int nu_, Function& fd_){
    updateVars(nx_, nu_, fd_);
};

void errorEstimator::updateVars(int nx_, int nu_, Function fd_){
    nx = nx_;
    nu = nu_;
    fd = fd_;
    px = std::vector<Polynomial>(nx);
    pu = std::vector<Polynomial>(nu);
};

double errorEstimator::getEstimate(double tf_, std::vector<double> tau, 
                                   std::vector<std::vector<double>>& xx,
                                   std::vector<std::vector<double>>& uu){
    int N = tau.size();
    tf = tf_;
    if (tau[0] > 0.0){
        tau.insert(tau.begin(), 0.);
    } else {
        N -= 1;
    }

    // Construct interpolant polynomial for later use
    constructInterpolant(tau, xx, uu);

    // Make sure xx and uu are equally long
    std::vector<double> uInit(nu);
    for (int k = 0; k < nu; k++){
        uInit[k] = pu[k](0.0);
    }
    uu.insert(uu.begin(), uInit);

    // Get new collocation points
    std::vector<double> tauHat = collocation_points(N+1);

    // Evaluate previous solution on new grid
    std::vector<std::vector<double>> xxHat(xx.size()+1, 
                                           std::vector<double>(nx));
    std::vector<std::vector<double>> uuHat(uu.size()+1, 
                                           std::vector<double>(nu));
    evaluateNewGrid(tau, tauHat, xx, xxHat, uu, uuHat);

    // Compute integrated values
    std::vector<std::vector<double>> yy(xxHat.size(), std::vector<double>(nx));
    getIntegratedValues(tauHat, xxHat, uuHat, yy);

    // Compute absolute errors and maximal state values
    std::vector<double> maxX(nx, 0.0);
    std::vector<std::vector<double>> err(xxHat.size(), 
                                         std::vector<double>(nx));
    for (int j = 0; j < xxHat.size(); j++){
        for (int k = 0; k < nx; k++){
            err[j][k] = abs(yy[j][k] - xxHat[j][k]);
            
            if (abs(xxHat[j][k]) > maxX[k]){ maxX[k] = abs(xxHat[j][k]);}
        }
    }

    // Compute relative errors
    double result = 0;
    for (int j = 0; j < xxHat.size(); j++){
        for (int k = 0; k < nx; k++){
            err[j][k] /= maxX[k];
            if (err[j][k] > result){ result = err[j][k];}
        }
    }
    
    return result;
};

void errorEstimator::evaluateNewGrid(std::vector<double>& tau, 
                                     std::vector<double>& tauHat,
                                     std::vector<std::vector<double>>& xx, 
                                     std::vector<std::vector<double>>& xxHat,
                                     std::vector<std::vector<double>>& uu, 
                                     std::vector<std::vector<double>>& uuHat){
    std::vector<double> etauHat = tauHat; etauHat.insert(etauHat.begin(), 0.);
    // for every state
    for (int k = 0; k < nx; k++){
        for (int j = 0; j < etauHat.size(); j++){
            xxHat[j][k] = px[k](etauHat[j]);
        }
    }

    // for every control
    for (int k = 0; k < nu; k++){
        for (int j = 0; j < etauHat.size(); j++){
            uuHat[j][k] = pu[k](etauHat[j]);
        }
    }
};

void errorEstimator::constructInterpolant(std::vector<double>& tau, 
                                        std::vector<std::vector<double>>& xx,
                                        std::vector<std::vector<double>>& uu){
    std::vector<double> etau = tau; 
    if (etau[0] > 0){etau.insert(etau.begin(), 0.);}
    Polynomial pLagrange;
    // for every state
    for (int k = 0; k < nx; k++){
        // Construct interplating polynomial by first constructing lagrange 
        // polynomials
        px[k] = 0;
        // for every lagrange polynomial
        for (int j = 0; j < etau.size(); j++){
            pLagrange = 1;
            // construct single Lagrange polynomial
            for (int r = 0; r < etau.size(); r++){
                if (r != j){
                    pLagrange *= Polynomial(-etau[r],1)/(etau[j]-etau[r]);
                }
            }
            // Add weighted lagrange polynomial to the interpolant
            pLagrange *= Polynomial(xx[j][k]);
            px[k] += pLagrange;
        }
        if (px[k].degree() < 0){
            px[k] = 0;
        }
    }

    // for every control
    for (int k = 0; k < nu; k++){
        // Construct interplating polynomial by first constructing lagrange 
        // polynomials
        pu[k] = 0;
        // for every lagrange polynomial
        for (int j = 1; j < etau.size(); j++){
            pLagrange = 1;
            // construct single Lagrange polynomial
            for (int r = 1; r < etau.size(); r++){
                if (r != j){
                    pLagrange *= Polynomial(-etau[r],1)/(etau[j]-etau[r]);
                }
            }
            // Add weighted lagrange polynomial to the interpolant
            pLagrange *= Polynomial(uu[j-1][k]);
            pu[k] += pLagrange;
        }
        if (pu[k].degree() < 0){
            pu[k] = 0;
        }
    }
};

// Function to construct an integration matrix
void errorEstimator::constructIntegrationMatrix(std::vector<double>& tauHat){
    // Get Integration matrix by constructing lagrange polynomials and 
    // integrating them
    std::vector<double> etauHat = tauHat; 
    if (etauHat[0] > 0){etauHat.insert(etauHat.begin(), 0.);}
    B = std::vector<std::vector<double>>(etauHat.size(), 
                                         std::vector<double>(etauHat.size()));
    Polynomial pLagrange;
    Polynomial ip;
    // for every lagrange polynomial
    for (int j = 0; j < etauHat.size(); j++){
        pLagrange = 1;
        
        // construct single Lagrange polynomial
        for (int r = 0; r < etauHat.size(); r++){
            if (r != j){
                pLagrange *= Polynomial(-etauHat[r],1)/(etauHat[j]-etauHat[r]);
            }
        }

        // Evaluate integration and store entries
        ip = pLagrange.anti_derivative();
        for (int i = 0; i < etauHat.size(); i++){
            B[i][j] = ip(etauHat[i]);
        }
    }
};

void errorEstimator::getIntegratedValues(std::vector<double>& tauHat, 
                                    std::vector<std::vector<double>>& xxHat,
                                    std::vector<std::vector<double>>& uuHat, 
                                    std::vector<std::vector<double>>& yy){

    // Construct the integration matrix
    constructIntegrationMatrix(tauHat);
    std::vector<double> etauHat = tauHat;
    if (etauHat[0] > 0){etauHat.insert(etauHat.begin(), 0.);}

    // Get dynamics function values
    std::vector<std::vector<double>> fvals(etauHat.size(), 
                                           std::vector<double>(nx));
    std::vector<const double*> arg(2);
    std::vector<double*> outputv(nx);
    std::vector<double> output(nx);
    for (int i = 0; i < nx; i++){ outputv[i] = &output[i];};
    
    // for every collocation point
    for (int j = 0; j < etauHat.size(); j++){
        
        // Evaluate dynamics
        arg[0] = &xxHat[j][0]; arg[1] = &uuHat[j][0];
        fd(arg, outputv);
        // store state derivatives
        for (int k = 0; k < nx; k++){
            // Here, a scaling is added to compensate the fact that we are 
            // evaluating on a grid of (0, 1] whereas the true interval has 
            // length tf.
            fvals[j][k] = output[k]*tf;
        }
    }

    // Compute integration using the integration matrix and the dynamics 
    // function values
    double integration;
    // for every collocation point

    for (int j = 0; j < etauHat.size(); j++){
        // for every state
        for (int k = 0; k < nx; k++){
            
            // Compute integration as a dot-product
            integration = 0.0;
            // for every lagrangian
            for (int i = 0; i < etauHat.size(); i++){
                integration += B[j][i]*fvals[i][k];
            }

            yy[j][k] = xxHat[0][k] + integration;
        }
    }
}
