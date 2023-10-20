
#include "../include/interfaceTester.hpp"
#include "../include/NLPInterface.hpp"
#include "IpTNLP.hpp"
#include <iostream>
#include <casadi/casadi.hpp>
#include <random>
#include <chrono>

using namespace casadi;
using namespace std;
using namespace chrono;


InterfaceTester::InterfaceTester(NLPInterface& interface_){
    interface = &interface_;
};

void InterfaceTester::get_starting_point(std::vector<double>& x){
    interface->get_starting_point(get_n(), true, &x[0], false, nullptr,
                                   nullptr, get_m(), false, nullptr);
    cout<<"initial guess = "<<endl;
    for (auto e : x){ cout<<e<<" ";} cout<<endl;
};

void InterfaceTester::eval_f(std::vector<double>& x){
    double res = 0;
    int n = get_n();
    interface->eval_f(n, &x[0], true, res);
    cout<<"objective function value = "<<res<<endl;
};

void InterfaceTester::eval_grad_f(std::vector<double>& x){
    std::vector<double> res(get_n());
    interface->eval_grad_f(get_n(), &x[0], true, &res[0]);
    cout<<"gradient value ="<<endl;
    for (auto elem : res){
        cout<<elem<<" ";
    }
    cout<<endl;
};

void InterfaceTester::eval_g(std::vector<double>& x){
    int n = get_n();
    int m = get_m();
    std::vector<double>res(m);
    std::vector<double>x_l(n);
    std::vector<double>x_u(n);
    std::vector<double>g_l(m);
    std::vector<double>g_u(m);
    interface->get_bounds_info(n, &x_l[0], &x_u[0], m, &g_l[0], &g_u[0]);

    interface->eval_g(n, &x[0], true, m, &res[0]);
    std::cout<<std::fixed<<std::setprecision(2);
    std::cout<<"constraint values ="<<endl;
    for (int i = 0; i < m; i++){
        if (g_l[i] >= 0){std::cout<<" "<<g_l[i];} else {std::cout<<g_l[i];}
        std::cout<<" <= ";
        if (res[i] >= 0){std::cout<<" "<<res[i];} else {std::cout<<res[i];}
        std::cout<<" <= ";
        if (g_u[i] >= 0){std::cout<<" "<<g_u[i];} else {std::cout<<g_u[i];}
        if (g_l[i] <= res[i] && res[i] <= g_u[i]){
            std::cout<<endl;
        } else {
            std::cout<<"   ("<<(g_l[i] <= res[i] && res[i] <= g_u[i])<<")"
                <<endl;
        }
    }
    std::cout<<std::defaultfloat;
};

void InterfaceTester::eval_jac_g(std::vector<double>& x){
    int n; int m; int nele_jac; int temp;
    Ipopt::TNLP::IndexStyleEnum index_style;
    interface->get_nlp_info(n, m, nele_jac, temp, index_style);

    std::vector<int> rows(nele_jac);
    std::vector<int> cols(nele_jac);
    std::vector<double> vals(nele_jac);

    interface->eval_jac_g(n, &x[0], true, m, nele_jac, &rows[0], &cols[0],
                          nullptr);
    interface->eval_jac_g(n, &x[0], true, m, nele_jac, nullptr, nullptr,
                          &vals[0]);
    
    // print the sparsity
    Sparsity sp = Sparsity(m, n);
    cout<<"m = "<<m<<", n = "<<n<<endl;
    for (int i = 0; i < nele_jac; i++){
        // cout<<"("<<rows[i]<<", "<<cols[i]<<")"<<endl;
        sp.add_nz(rows[i], cols[i]);
    }

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            if (sp.has_nz(i,j)){
                cout<<"X ";
            } else {
                cout<<". ";
            }
        }
        cout<<endl;
    }
    cout<<endl;

    // print the values
    // for (int i = 0; i < nele_jac; i++){
    //     cout<<"["<<i<<"] ("<<rows[i]<<", "<<cols[i]<<") "<<vals[i]<<endl;
    // }
};

void InterfaceTester::eval_h(std::vector<double>& x,
                             std::vector<double>& lambda){
    int n; int m; int nele_hess; int temp;
    Ipopt::TNLP::IndexStyleEnum index_style;
    interface->get_nlp_info(n, m, temp, nele_hess, index_style);

    std::vector<int> rows(nele_hess);
    std::vector<int> cols(nele_hess);
    std::vector<double> vals(nele_hess);

    interface->eval_h(n, &x[0], true, 1.0, m, &lambda[0], true, nele_hess,
                      &rows[0], &cols[0], nullptr);

    interface->eval_h(n, &x[0], true, 1.0, m, &lambda[0], true, nele_hess,
                      nullptr, nullptr, &vals[0]);
    
    // print the sparsity
    Sparsity sp = Sparsity(n, n);
    for (int i = 0; i < nele_hess; i++){
        sp.add_nz(rows[i], cols[i]);
    }

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (sp.has_nz(i,j)){
                cout<<"X ";
            } else {
                cout<<". ";
            }
        }
        cout<<endl;
    }
    cout<<endl;

    // print the values
    for (int i = 0; i < nele_hess; i++){
        cout<<"("<<rows[i]<<", "<<cols[i]<<") "<<vals[i]<<endl;
    }
};

std::vector<double> InterfaceTester::computeFunctionEvaluationTimes(){
    // query interface information
    int n; int m; int nele_jac; int nele_hess;
    Ipopt::TNLP::IndexStyleEnum index_style;
    interface->get_nlp_info(n, m, nele_jac, nele_hess, index_style);

    // define input vectors
    int nb_calls = 1000;
    std::vector<std::vector<double>> x(nb_calls, std::vector<double>(n));
    std::vector<std::vector<double>> lambda(nb_calls, std::vector<double>(m));

    // construct output containers    
    double obj;
    std::vector<double> grad(n);
    std::vector<double> g(m);
    std::vector<double> jac_vals(nele_jac);
    std::vector<double> hess_vals(nele_hess);
    std::vector<double> computation_times(5);

    // construct randomizer and random values
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-10.0, 10.0);
    for (int i = 0; i < nb_calls; i++){
        for (int j = 0; j < n; j++){
            x[i][j] = dis(gen);
        }
        for (int j = 0; j < m; j++){
            lambda[i][j] = dis(gen);
        }
    }

    cout<<"starting f"<<endl;
    auto start = high_resolution_clock::now();
    for (int i = 0; i < nb_calls; i++){
        interface->eval_f(n, &x[i][0], true, obj);
    }
    auto stop = high_resolution_clock::now();
    computation_times[0] = double(duration_cast<microseconds>(stop-start).count())/1000.0/nb_calls;    

    cout<<"starting grad_f"<<endl;
    start = high_resolution_clock::now();
    for (int i = 0; i < nb_calls; i++){
        interface->eval_grad_f(n, &x[i][0], true, &grad[0]);
    }
    stop = high_resolution_clock::now();
    computation_times[1] = double(duration_cast<microseconds>(stop-start).count())/1000.0/nb_calls;    

    cout<<"starting g"<<endl;
    start = high_resolution_clock::now();
    for (int i = 0; i < nb_calls; i++){
        interface->eval_g(n, &x[i][0], true, m, &g[0]);
    }
    stop = high_resolution_clock::now();
    computation_times[2] = double(duration_cast<microseconds>(stop-start).count())/1000.0/nb_calls;    

    cout<<"starting jac"<<endl;
    start = high_resolution_clock::now();
    for (int i = 0; i < nb_calls; i++){
        interface->eval_jac_g(n, &x[i][0], true, m, nele_jac, nullptr, nullptr,
                          &jac_vals[0]);
    }
    stop = high_resolution_clock::now();
    computation_times[3] = double(duration_cast<microseconds>(stop-start).count())/1000.0/nb_calls;    

    cout<<"starting hess"<<endl;
    start = high_resolution_clock::now();
    for (int i = 0; i < nb_calls; i++){
        interface->eval_h(n, &x[i][0], true, 1.0, m, &lambda[i][0], true, nele_hess,
                          nullptr, nullptr, &hess_vals[0]);
    }
    stop = high_resolution_clock::now();
    computation_times[4] = double(duration_cast<microseconds>(stop-start).count())/1000.0/nb_calls;    

    for (auto e : computation_times){cout<<e<<" ";} cout<<endl;
    return computation_times;
};


int InterfaceTester::get_n(){
    int n; int m; int nnz_jac_g; int nnz_h_lag;
    Ipopt::TNLP::IndexStyleEnum index_style;
    interface->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, 
                            index_style);
    return n;
};

int InterfaceTester::get_m(){
    int n; int m; int nnz_jac_g; int nnz_h_lag;
    Ipopt::TNLP::IndexStyleEnum index_style;
    interface->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, 
                            index_style);
    return m;
};
