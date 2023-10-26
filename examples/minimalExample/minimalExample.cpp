#include "core.hpp"
#include "makeBuildingBlocks.hpp"
#include "plotter.hpp"
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;


// Function to generate a .csv file with the solutions that will be put in 
// the folder plotting_data
void writeSolutionToFile(AdaptiveNLP& nlp, std::vector<double>& sol, 
                         int file_nb){
    int nx = nlp.getNx();
    int nu = nlp.getNu();
    int N = nlp.getN();
    int final_ind = nlp.getFinalInd();
    
    double t_sol = sol[0];

    std::string file_path = __FILE__;
    std::string dir_path = file_path.substr(0, file_path.rfind("/"));
    std::ofstream file(dir_path + "/plotting_data/solution_" + 
                       std::to_string(file_nb) + ".csv");
    if (file.is_open()){
        file << "t";
        for (int i = 0; i < nx; i++){
            file << ",x_" + std::to_string(i);
        }
        for (int i = 0; i < nu; i++){
            file << ",u_" + std::to_string(i);
        }
        file << "\n";

        // loop over all time-steps
        int k = 0;
        while (k != final_ind){
            // write time-stamp
            file << std::to_string(t_sol * nlp.getTimeFromIndex(k));
            
            // write state and control values
            for (int i = 0; i < nx + nu; i++){
                file << "," <<
                        std::to_string(sol[nlp.getFreeTime() + nx*k + 
                                           nu*k - nu*(k > final_ind) + i]);
            }
            file << "\n";
            k = nlp.getNext(k);
        }
        // write final state values
        file << std::to_string(t_sol * nlp.getTimeFromIndex(final_ind));
        for (int i = 0; i < nx; i++){
            file << "," <<
                std::to_string(sol[nlp.getFreeTime() + nx*final_ind + 
                                    nu*final_ind - nu*(k > final_ind) + i]);
        }
        // put placeholder control values 
        for (int i = 0; i < nu; i++){
            file << "," << std::to_string(-1);
        }

        file.close();
    }
}

int main(){
    // Construct building blocks
    BuildingBlocks blocks = minimalBlocks();

    // Create AdaptiveNLP instance
    AdaptiveNLP adaptiveNLP = AdaptiveNLP(blocks, 1.0, 60, 5);

    // Initialize time-steps
    int nb_intervals = 20;
    std::vector<int> nks(nb_intervals, 2);
    adaptiveNLP.initTimeSteps(nks);

    // Solve the NLP
    adaptiveNLP.solveNlp({{"p_g0", std::vector<double>{0.0, 0.0, 0.0, 0.0}}, 
                         {"p_gT", std::vector<double>{3.0, 1.0, 0.0, 0.0}}});
    std::vector<double> sol = adaptiveNLP.getSolution();
    writeSolutionToFile(adaptiveNLP, sol, 0);

    // Add no-collision constraint to all points
    std::vector<int> k_to_add(adaptiveNLP.getN());
    for (int i = 0; i <= adaptiveNLP.getN(); i++){
        if (i != adaptiveNLP.getFinalInd()){
            k_to_add[i - (i > adaptiveNLP.getFinalInd())] = i;
        }
    }

    // Add an obstacle at position (1.5, 0.4) with radius 0.3
    adaptiveNLP.addExtraConstraint(k_to_add, {0}, {{1.5, 0.4, 0.3}});

    // Let's solve again, this time providing our previous solution as an 
    // intitial guess
    adaptiveNLP.solveNlp({{"p_g0", std::vector<double>{0.0, 0.0, 0.0, 0.0}}, 
                         {"p_gT", std::vector<double>{3.0, 1.0, 0.0, 0.0}}}, 
                         sol);
    sol = adaptiveNLP.getSolution();
    writeSolutionToFile(adaptiveNLP, sol, 1);

    // Now, let's refine the time-grid in the first interval
    // If a uniform refinement is fine, we can leave out the third argument 
    // which is a vector of time-stamps for the inserted time-steps
    int N_old = adaptiveNLP.getN();
    adaptiveNLP.changeIntervalDiscretization(9, {2, 2, 2});
    adaptiveNLP.changeIntervalDiscretization(10, {2, 2, 2});
    adaptiveNLP.changeIntervalDiscretization(11, {2, 2, 2});

    // We still have to add the extra constraint to these new time-steps
    std::vector<int> k_to_add_new(adaptiveNLP.getN() - N_old);
    for (int i = 0; i < k_to_add_new.size(); i++){
        k_to_add_new[i] = adaptiveNLP.getFinalInd() + 1 +i;
    }
    adaptiveNLP.addExtraConstraint(k_to_add_new, {0}, {{1.5, 0.4, 0.3}});

    // Let's solve the NLP again
    try{
        adaptiveNLP.solveNlp({{"p_g0", std::vector<double>{0.0, 0.0, 0.0, 0.0}}, 
                            {"p_gT", std::vector<double>{3.0, 1.0, 0.0, 0.0}}}, 
                            sol);
        sol = adaptiveNLP.getSolution();
        writeSolutionToFile(adaptiveNLP, sol, 2);
    } catch (const std::exception& exc) {
        std::cout<<"Exception caught: "<<exc.what()<<endl; 
        // AdaptiveNLP throws an exception because the provided initial guess
        // is not of correct length due to the refinements we've made. 
        // Let's just omit the initial guess. In that case, AdaptiveNLP will
        // consider the building blocks provided to evaluate an initial guess.
        adaptiveNLP.solveNlp({{"p_g0", std::vector<double>{0.0, 0.0, 0.0, 0.0}}, 
                            {"p_gT", std::vector<double>{3.0, 1.0, 0.0, 0.0}}});
        sol = adaptiveNLP.getSolution();
        writeSolutionToFile(adaptiveNLP, sol, 2);
    }

    return 0;
}