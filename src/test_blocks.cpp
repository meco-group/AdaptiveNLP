#include "../include/makeBuildingBlocks.hpp"
#include "../include/bookkeeper.hpp"
#include "../include/NLPInterface.hpp"
#include "../include/adaptiveNLP.hpp"
#include "../include/plotter.hpp"
#include <iostream>
#include <vector>

using namespace std;

int main(){
    std::vector<double> corridor_1;
    std::vector<double> corridor_2;
    BuildingBlocks blocks = makeBlocksCorridors(3, corridor_1, corridor_2);
    // BuildingBlocks blocks = makeSimpleBlocks();
    // BuildingBlocks blocks = makeSimpleBlocksFreeTime();
    cout<<"Building blocks defined"<<endl;

    Plotter plotter = Plotter(blocks.getFreeTime(), blocks.get_nx(), 
                              blocks.get_nu());

    AdaptiveNLP adaptiveNLP = AdaptiveNLP(blocks, 5.0, 60, 5, &plotter);
    cout<<"adaptiveNLP is constructed"<<endl;

    int nb_intervals = 15;
    std::vector<int> nks(nb_intervals, 3);
    adaptiveNLP.initTimeSteps(nks);
    cout<<"initialized time steps"<<endl;

    adaptiveNLP.solveNlp({{"p_g0", std::vector<double>{0.0, 0.0, 0.0, 0.0}}, 
                         {"p_gT", std::vector<double>{2.0, 2.1, -0.2, 0.2}}});
    std::vector<double> sol = adaptiveNLP.getSolution();

    std::vector<int> all(nb_intervals*2);
    for (int i = 0; i < all.size(); i++){all[i] = i;}
    std::vector<int> constraint_inds(all.size(), 1);
    std::vector<std::vector<double>> params(all.size(), std::vector<double>(3));
    for (int i = 0; i < all.size(); i++){params[i] = {1.0, 0.0, 0.5};};
    adaptiveNLP.addExtraConstraint(all, constraint_inds, params);
    plotter.addCircle(1.0, 0.0, 0.5);

    std::vector<double> corridor = {2.0, 2.1, -0.2, 0.2};
    plotter.addCorridor(corridor);

    adaptiveNLP.solveNlp({{"p_g0", std::vector<double>{0.0, 0.0, 0.0, 0.0}}, 
                         {"p_gT", std::vector<double>{2.0, 2.1, -0.2, 0.2}}}, 
                         sol);

    // cout<<"appending..."<<endl;
    adaptiveNLP.appendTimeSteps({3, 3});
    // adaptiveNLP.appendTimeSteps({3});
    // cout<<"done!"<<endl;

    adaptiveNLP.addTimeSteps({3});
    // adaptiveNLP.addTimeSteps({3, 3, 3, 3, 3, 3, 3});

    adaptiveNLP.solveNlp({{"p_g0", std::vector<double>{0.0, 0.0, 0.0, 0.0}}, 
                         {"p_gT", std::vector<double>{2.0, 2.1, -0.2, 0.2}}});
    // adaptiveNLP.solveNlp({});

    // adaptiveNLP.testInterface();


    return 0;
}