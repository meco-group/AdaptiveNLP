#ifndef __INTERFACE_TESTER__
#define __INTERFACE_TESTER__
#include "NLPInterface.hpp"
#include <iostream>
#include "IpTNLP.hpp"

using namespace std;

class InterfaceTester{
    public:
        InterfaceTester(){};

        InterfaceTester(NLPInterface& interface_);

        void get_starting_point(std::vector<double>& x);

        void eval_f(std::vector<double>& x);

        void eval_grad_f(std::vector<double>& x);

        void eval_g(std::vector<double>& x);

        void eval_jac_g(std::vector<double>& x);

        void eval_h(std::vector<double>& x, std::vector<double>& lambda);

        std::vector<double> computeFunctionEvaluationTimes();

    private:
        int get_n();

        int get_m();

        NLPInterface* interface;


};

#endif