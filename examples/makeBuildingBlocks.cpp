#include "core.hpp"
#include <casadi/casadi.hpp>
#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace casadi;
using namespace std;

BuildingBlocks minimalBlocks(){
    int nx = 4;
    int nu = 2;

    // Helper variables
    SX xk = SX::sym("xk", nx, 1);
    SX uk = SX::sym("uk", nu, 1);
    SX t = SX::sym("t", 1, 1);
    SX p = SX::sym("p", 0, 0);
    std::vector<SX> arg = {xk, t, p};

    BuildingBlocks myBlocks = BuildingBlocks();
    myBlocks.add_sizes(nx, nu, true);

    // Initial cost
    Function Phi_0("Phi_0", arg, {t});
    Function Phi_0_J("Phi_0_J", arg, {jacobian(Phi_0(arg)[0], vertcat(t, xk))});
    Function Phi_0_H("Phi_0_H", arg, {hessian(Phi_0(arg)[0], vertcat(t, xk))});
    myBlocks.add_Phi_0(Phi_0, Phi_0_J, Phi_0_H);

    // Final cost
    Function Phi_f("Phi_f", arg, {0.0});
    Function Phi_f_J("Phi_f_J", arg, {jacobian(Phi_f(arg)[0], vertcat(t, xk))});
    Function Phi_f_H("Phi_f_H", arg, {hessian(Phi_f(arg)[0], vertcat(t, xk))});
    myBlocks.add_Phi_f(Phi_f, Phi_f_J, Phi_f_H);

    // Stage cost
    SX dt = SX::sym("dt", 1, 1);
    std::vector<SX> arg_phi = {xk, uk, dt, t, p};
    Function phi("phi", arg_phi, {0.0});
    Function phi_J("phi_J", arg_phi, {jacobian(phi(arg_phi)[0], 
                                               vertcat(t, xk, uk))});
    Function phi_H("phi_H", arg_phi, {hessian(phi(arg_phi)[0], 
                                              vertcat(t, xk, uk))});
    myBlocks.add_phi(phi, phi_J, phi_H);

    // Initial constraint
    SX lambda_0 = SX::sym("lambda_0", 5, 1);
    SX p_g0 = SX::sym("p_g0", 4, 1);
    std::vector<SX> arg_g0 = {xk, t, p_g0};
    Function g0("g0", arg_g0, {vertcat(xk-p_g0, t)});
    std::vector<double> g0_lb = {0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> g0_ub = {0.0, 0.0, 0.0, 0.0, 20.0};
    Function g0_J("g0_J", arg_g0, {jacobian(g0(arg_g0)[0], vertcat(t, xk))});
    Function g0_H("g0_H", {xk, t, lambda_0, p_g0},
                  {hessian(dot(lambda_0, g0(arg_g0)[0]), vertcat(t, xk))});
    myBlocks.add_g0(g0, g0_J, g0_H, g0_lb, g0_ub, 4);

    // Final constraint
    SX lambda_T = SX::sym("lambda_T", 3, 1);
    SX p_gT = SX::sym("p_gT", 3, 1);
    std::vector<SX> arg_gT = {xk, p_gT};
    Function gT("gT", arg_gT, {xk(Slice(0,3))-p_gT});
    std::vector<double> gT_lb = {0.0, 0.0, 0.0};
    std::vector<double> gT_ub = {0.0, 0.0, 0.0};
    Function gT_J("gT_J", arg_gT, {jacobian(gT(arg_gT)[0], xk)});
    Function gT_H("gT_H", {xk, lambda_T, p_gT},
                             {hessian(dot(lambda_T, gT(arg_gT)[0]), xk)});
    myBlocks.add_gT(gT, gT_J, gT_H, gT_lb, gT_ub, 4);

    // Fixed constraints
    std::vector<SX> arg_g_fixed = {xk, uk, p};
    Function g_fixed("g_fixed", arg_g_fixed, {vertcat(xk(2), xk(3), uk)});
    std::vector<double> g_fixed_lb = {0.0, -M_PI/2.0, -1.0, -M_PI/3.0};
    std::vector<double> g_fixed_ub = {2.0, M_PI/2.0, 1.0, M_PI/3.0};
    int nb_g_fixed_ = g_fixed.size1_out(0);
    SX lambda_g_fixed = SX::sym("lambda_g_fixed", nb_g_fixed_, 1);
    Function g_fixed_J("g_fixed_J", arg_g_fixed,
                                  {jacobian(g_fixed(arg_g_fixed)[0],
                                            vertcat(xk, uk))});
    Function g_fixed_H("g_fixed_H", {xk, uk, lambda_g_fixed, p},
                                  {hessian(dot(lambda_g_fixed,
                                               g_fixed(arg_g_fixed)[0]),
                                           vertcat(xk, uk))});
    myBlocks.add_g_fixed(g_fixed, g_fixed_J, g_fixed_H, g_fixed_lb,
                         g_fixed_ub);

    // Flexible constraints
    //  obstacle-avoidance constraint
    SX p_extra = SX::sym("p_extra", 3, 1);
    std::vector<SX> arg_g_extra = {xk, uk, p_extra};
    Function g_extra("g_extra", arg_g_extra, {pow(xk(0)-p_extra(0), 2) +
                                              pow(xk(1)-p_extra(1), 2) -
                                              pow(p_extra(2), 2)});
    std::vector<double> lb_extra = {0.0};
    std::vector<double> ub_extra = {1.0e10};
    SX lambda_g_extra = SX::sym("lambda_extra", g_extra.size1_out(0), 1);
    Function g_extra_J("g_extra_J", arg_g_extra,
                                    {jacobian(g_extra(arg_g_extra)[0],
                                              vertcat(xk, uk))});
    Function g_extra_H("g_extra_H", 
                                  {xk, uk, lambda_g_extra, p_extra},
                                  {hessian(dot(lambda_g_extra,
                                               g_extra(arg_g_extra)[0]),
                                            vertcat(xk, uk))});

    myBlocks.add_g_extra({g_extra}, {g_extra_J}, {g_extra_H}, {lb_extra}, 
                         {ub_extra}, {3});

    // Discretization constraints
    SX lambda_disc2 = SX::sym("lambda_disc", 4, 1);
    SX dt2 = SX::sym("dt", 1, 1);
    SX xx2 = SX::sym("xx", nx, 2);
    SX uu2 = SX::sym("uu", nu, 1);
    SX rhs2 = vertcat(vertcat(xx2(2,0)*cos(xx2(3,0)), xx2(2,0)*sin(xx2(3,0))),
                        vertcat(uu2(0,0), xx2(2.0)*tan(uu2(1,0))));
    const std::vector<SX> arg_disc2 = {xx2, uu2, dt2, t, p};
    Function g_disc2("g_disc", arg_disc2, 
                     {xx2(Slice(),1) - (xx2(Slice(),0) + t*dt2*rhs2)});
    Function g_disc_J2("g_disc_J", arg_disc2,
                       {jacobian(g_disc2(arg_disc2)[0],
                                 vertcat(vertcat(t, xx2(Slice(),0)),
                                         uu2(Slice(),0),
                                         xx2(Slice(),1)))});
    Function g_disc_H2("g_disc_H", {xx2, uu2, dt2, t, lambda_disc2, p},
                       {hessian(dot(lambda_disc2, g_disc2(arg_disc2)[0]),
                                    vertcat(vertcat(t, xx2(Slice(), 0)),
                                            uu2(Slice(),0),
                                            xx2(Slice(), 1)))});

    SX lambda_disc3 = SX::sym("lambda_disc", 8, 1);
    SX dt3 = SX::sym("dt", 2, 1);
    SX xx3 = SX::sym("xx", nx, 3);
    SX uu3 = SX::sym("uu", nu, 2);
    SX rhs31 = vertcat(vertcat(xx3(2,0)*cos(xx3(3,0)),
                                xx3(2,0)*sin(xx3(3,0))),
                        vertcat(uu3(0,0), xx3(2,0)*tan(uu3(1,0))));
    SX rhs32 = vertcat(vertcat(xx3(2,1)*cos(xx3(3,1)),
                                xx3(2,1)*sin(xx3(3,1))),
                        vertcat(uu3(0,1), xx3(2,1)*tan(uu3(1,1))));
    const std::vector<SX> arg_g_disc3 = {xx3, uu3, dt3, t, p};
    Function g_disc3("g_disc", arg_g_disc3,
                     {vertcat(xx3(Slice(),1) -
                              (xx3(Slice(),0) + t*dt3(0)*rhs31),
                              xx3(Slice(),2) -
                              (xx3(Slice(),1) + t*dt3(1)*rhs32))});
    Function g_disc_J3("g_disc_J", arg_g_disc3,
                       {jacobian(g_disc3(arg_g_disc3)[0], 
                                 vertcat(t, vertcat(xx3(Slice(),0),
                                                 uu3(Slice(),0),
                                                 xx3(Slice(),1)),
                                         vertcat(uu3(Slice(),1),
                                                 xx3(Slice(),2))))});
    Function g_disc_H3("g_disc_H", {xx3, uu3, dt3, t, lambda_disc3, p},
                        {hessian(dot(lambda_disc3, g_disc3(arg_g_disc3)[0]),
                                    vertcat(t, vertcat(xx3(Slice(),0),
                                                    uu3(Slice(),0),
                                                    xx3(Slice(),1)),
                                            vertcat(uu3(Slice(),1),
                                                    xx3(Slice(),2))))});
    myBlocks.add_g_disc({g_disc2, g_disc3}, {g_disc_J2, g_disc_J3}, 
                        {g_disc_H2, g_disc_H3});

    // Initialization
    double T = 5.0;
    Function x_init("x_init", {t}, {vertcat(vertcat(3*t, t),
                                                       SX::zeros(2,1))});
    Function u_init("u_init", {t}, {SX::zeros(2, 1)});
    myBlocks.add_inits(x_init, u_init);

    return myBlocks;
};

BuildingBlocks makeBlocksCorridors(int nb_coll,
                                   std::vector<double>& corridor_1,
                                   std::vector<double>& corridor_2){
    int nx = 4;
    int nu = 2;
    corridor_1 = {-10.0, 5.0, 0.5, 4.0};
    corridor_2 = {-2.0, 35.0, 2.0, 6.0};
    double alpha = 1.0;

    // Helper variables
    SX xk = SX::sym("xk", nx, 1);
    SX uk = SX::sym("uk", nu, 1);
    SX p = SX::sym("p", 0, 0);
    std::vector<SX> arg = {xk, p};

    BuildingBlocks myBlocks = BuildingBlocks();
    myBlocks.add_sizes(nx, nu);

    // Initial cost
    Function Phi_0("Phi_0", {xk, p}, {0.0});
    Function Phi_0_J("Phi_0_J", {xk, p}, {jacobian(Phi_0(arg)[0], xk)});
    Function Phi_0_H("Phi_0_H", {xk, p}, {hessian(Phi_0(arg)[0], xk)});
    myBlocks.add_Phi_0(Phi_0, Phi_0_J, Phi_0_H);

    // Final cost
    Function Phi_f("Phi_f", {xk, p}, {-alpha * xk(0)});
    Function Phi_f_J("Phi_f_J", {xk, p}, {jacobian(Phi_f(arg)[0], xk)});
    Function Phi_f_H("Phi_f_H", {xk, p}, {hessian(Phi_f(arg)[0], xk)});
    myBlocks.add_Phi_f(Phi_f, Phi_f_J, Phi_f_H);

    // Stage cost
    SX dt = SX::sym("dt", 1, 1);
    std::vector<SX> arg_phi = {xk, uk, dt, p};
    Function phi("phi", {xk, uk, dt, p}, {dt * (uk(0)*uk(0) + uk(1)*uk(1))});
    Function phi_J("phi_J", {xk, uk, dt, p}, 
                   {jacobian(phi(arg_phi)[0], vertcat(xk, uk))});
    Function phi_H("phi_H", {xk, uk, dt, p}, 
                    {hessian(phi(arg_phi)[0], vertcat(xk, uk))});
    myBlocks.add_phi(phi, phi_J, phi_H);

    // Initial constraint
    SX lambda_0 = SX::sym("lambda_0", 4, 1);
    SX p_g0 = SX::sym("p_g0", 4, 1);
    SX init_cond = SX::zeros(4, 1);
    init_cond(0) = xk(0) - p_g0(0); init_cond(1) = xk(1) - p_g0(1);
    init_cond(2) = xk(2) - p_g0(2); init_cond(3) = xk(3) - p_g0(3);
    std::vector<SX> arg_g0 = {xk, p_g0};

    Function g0("g0", {xk, p_g0}, {init_cond});
    std::vector<double> g0_lb = {0.0, 0.0, 0.0, 0.0};
    std::vector<double> g0_ub = {0.0, 0.0, 0.0, 0.0};
    Function g0_J("g0_J", {xk, p_g0}, {jacobian(g0(arg_g0)[0], xk)});
    Function g0_H("g0_H", {xk, lambda_0, p_g0},
                 {hessian(dot(lambda_0, g0(arg_g0)[0]), xk)});
    myBlocks.add_g0(g0, g0_J, g0_H, g0_lb, g0_ub, 4);

    // Final constraint
    SX lambda_T = SX::sym("lambda_T", 4, 1);
    SX p_gT = SX::sym("p_gT", 4, 1);
    std::vector<SX> arg_gT = {xk, p_gT};
    Function gT("gT", {xk, p_gT}, {vertcat(vertcat(xk(0)-p_gT(0),
                                                   p_gT(1)-xk(0)),
                                           vertcat(xk(1)-p_gT(2),
                                                   p_gT(3)-xk(1)))});
    std::vector<double> gT_lb = {0.0, 0.0, 0.0, 0.0};
    std::vector<double> gT_ub = {1.0e10, 1.0e10, 1.0e10, 1.0e10};
    Function gT_J("gT_J", {xk, p_gT}, {jacobian(gT(arg_gT)[0], xk)});
    Function gT_H("gT_H", {xk, lambda_T, p_gT},
                  {hessian(dot(lambda_T, gT(arg_gT)[0]), xk)});
    myBlocks.add_gT(gT, gT_J, gT_H, gT_lb, gT_ub, 4);

    // Fixed constraints
    Function g_fixed("g_fixed", {xk, uk, p}, {vertcat(xk(2), xk(3), uk)});
    std::vector<double> g_fixed_lb = {0.0, -M_PI/2.0, -1.0, -M_PI/3.0};
    std::vector<double> g_fixed_ub = {2.0, M_PI/2.0, 1.0, M_PI/3.0};
    int nb_g_fixed_ = g_fixed.size1_out(0);
    SX lambda_g_fixed = SX::sym("lambda_g_fixed", nb_g_fixed_, 1);
    std::vector<SX> arg_g_fixed = {xk, uk, p};
    Function g_fixed_J("g_fixed_J", {xk, uk, p}, 
                       {jacobian(g_fixed(arg_g_fixed)[0], vertcat(xk, uk))});
    Function g_fixed_H("g_fixed_H", {xk, uk, lambda_g_fixed, p},
                       {hessian(dot(lambda_g_fixed, g_fixed(arg_g_fixed)[0]),
                                vertcat(xk, uk))});
    myBlocks.add_g_fixed(g_fixed, g_fixed_J, g_fixed_H, g_fixed_lb,
                         g_fixed_ub);


    // Flexible constraints
    //  1: corridor constraint
    SX p_extra_1 = SX::sym("p_extra_1", 4, 1);
    Function g_extra_1("g_extra_1", {xk, uk, p_extra_1},
                                  {vertcat(vertcat(xk(0) - p_extra_1(0),
                                                   p_extra_1(1) - xk(0)),
                                           vertcat(xk(1) - p_extra_1(2),
                                                   p_extra_1(3) - xk(1)))});
    std::vector<double> lb_1 = {0.0, 0.0, 0.0, 0.0};
    std::vector<double> ub_1 = {1.0e10, 1.0e10, 1.0e10, 1.0e10};
    SX lambda_g_extra_1 = SX::sym("lambda_extra_1", g_extra_1.size1_out(0), 1);
    std::vector<SX> arg_g_extra_1 = {xk, uk, p_extra_1};
    Function g_extra_J_1("g_extra_J_1", {xk, uk, p_extra_1},
                                    {jacobian(g_extra_1(arg_g_extra_1)[0],
                                              vertcat(xk, uk))});
    Function g_extra_H_1("g_extra_H_1", {xk, uk, p_extra_1},
                                    {hessian(
                                        dot(lambda_g_extra_1,
                                            g_extra_1(arg_g_extra_1)[0]),
                                        vertcat(xk, uk))});

    //  2: obstacle-avoidance constraint
    SX p_extra_2 = SX::sym("p_extra_2", 3, 1);
    Function g_extra_2("g_extra_2", {xk, uk, p_extra_2},
                                  {pow(xk(0)-p_extra_2(0), 2) +
                                   pow(xk(1)-p_extra_2(1), 2) -
                                   pow(p_extra_2(2), 2)});
    std::vector<double> lb_2 = {0.0};
    std::vector<double> ub_2 = {1.0e10};
    SX lambda_g_extra_2 = SX::sym("lambda_extra_2", g_extra_2.size1_out(0), 1);
    std::vector<SX> arg_g_extra_2 = {xk, uk, p_extra_2};
    Function g_extra_J_2("g_extra_J_2", {xk, uk, p_extra_2},
                         {jacobian(g_extra_2(arg_g_extra_2)[0], 
                                   vertcat(xk, uk))});
    Function g_extra_H_2("g_extra_H_2", {xk, uk, lambda_g_extra_2, p_extra_2},
                         {hessian(dot(lambda_g_extra_2, 
                                      g_extra_2(arg_g_extra_2)[0]),
                                  vertcat(xk, uk))});

    // 3: safety constraint
    SX p_extra_3 = SX::sym("p_extra_3", 4, 1);
    SX R = p_extra_3(1); 
    double vmax = 2.0;
    SX d = p_extra_3(0);
    SX a = (vmax - d)/pow(R, 3);
    SX b = -(vmax - d + 2*a*pow(R, 3))/pow(R, 2);
    SX c = -(3*a*pow(R, 2) + 2*b*R);
    SX r = sqrt(pow(xk(0)-p_extra_3(2), 2) +
           pow(xk(1)-p_extra_3(3), 2));
    Function g_extra_3("g_extra_3", {xk, uk, p_extra_3},
                                  {a*pow(r, 3) + b*pow(r, 2) +
                                  c*r + d - xk(2)});
    std::vector<double> lb_3 = {0.0};
    std::vector<double> ub_3 = {1.0e10};
    SX lambda_g_extra_3 = SX::sym("lambda_extra_3", g_extra_3.size1_out(0), 1);
    std::vector<SX> arg_g_extra_3 = {xk, uk, p_extra_3};
    Function g_extra_J_3("g_extra_J_3", {xk, uk, p_extra_3},
                         {jacobian(g_extra_3(arg_g_extra_3)[0],
                                   vertcat(xk, uk))});
    Function g_extra_H_3("g_extra_H_3", {xk, uk, lambda_g_extra_3, p_extra_3},
                         {hessian(dot(lambda_g_extra_3,
                                      g_extra_3(arg_g_extra_3)[0]),
                                  vertcat(xk, uk))});

    myBlocks.add_g_extra({g_extra_1, g_extra_2, g_extra_3},
                         {g_extra_J_1, g_extra_J_2, g_extra_J_3},
                         {g_extra_H_1, g_extra_H_2, g_extra_H_3},
                         {lb_1, lb_2, lb_3}, {ub_1, ub_2, ub_3},
                         {4, 3, 4});

    // Discretization constraints
    SX lambda_disc2 = SX::sym("lambda_disc", 4, 1);
    SX dt2 = SX::sym("dt", 1, 1);
    SX xx2 = SX::sym("xx", nx, 2);
    SX uu2 = SX::sym("uu", nu, 1);
    SX rhs2 = vertcat(vertcat(xx2(2,0)*cos(xx2(3,0)), xx2(2,0)*sin(xx2(3,0))),
                      vertcat(uu2(0,0), xx2(2.0)*tan(uu2(1,0))));
    const std::vector<SX> arg_disc2 = {xx2, uu2, dt2, p};
    Function g_disc2("g_disc", {xx2, uu2, dt2, p},
                                {xx2(Slice(),1) - (xx2(Slice(),0) + dt2*rhs2)});
    Function g_disc_J2("g_disc_J", {xx2, uu2, dt2, p},\
                                    {jacobian(g_disc2(arg_disc2)[0],
                                            vertcat(xx2(Slice(),0),
                                                    uu2(Slice(),0),
                                                    xx2(Slice(),1)))});
    Function g_disc_H2("g_disc_H", {xx2, uu2, dt2, lambda_disc2, p},\
                                    {hessian(dot(lambda_disc2,
                                                g_disc2(arg_disc2)[0]),
                                            vertcat(xx2(Slice(), 0),
                                                    uu2(Slice(),0),
                                                    xx2(Slice(), 1)))});

    SX lambda_disc3 = SX::sym("lambda_disc", 8, 1);
    SX dt3 = SX::sym("dt", 2, 1);
    SX xx3 = SX::sym("xx", nx, 3);
    SX uu3 = SX::sym("uu", nu, 2);
    SX rhs31 = vertcat(vertcat(xx3(2,0)*cos(xx3(3,0)),
                                xx3(2,0)*sin(xx3(3,0))),
                        vertcat(uu3(0,0), xx3(2,0)*tan(uu3(1,0))));
    SX rhs32 = vertcat(vertcat(xx3(2,1)*cos(xx3(3,1)),
                                xx3(2,1)*sin(xx3(3,1))),
                        vertcat(uu3(0,1), xx3(2,1)*tan(uu3(1,1))));
    const std::vector<SX> arg_g_disc3 = {xx3, uu3, dt3, p};
    Function g_disc3("g_disc", {xx3, uu3, dt3, p},
                                {vertcat(xx3(Slice(),1) -
                                        (xx3(Slice(),0) + dt3(0)*rhs31),
                                        xx3(Slice(),2) -
                                        (xx3(Slice(),1) + dt3(1)*rhs32))});
    Function g_disc_J3("g_disc_J", {xx3, uu3, dt3, p},
                                    {jacobian(g_disc3(arg_g_disc3)[0], 
                                            vertcat(vertcat(xx3(Slice(),0),
                                                            uu3(Slice(),0),
                                                            xx3(Slice(),1)),
                                                    vertcat(uu3(Slice(),1),
                                                            xx3(Slice(),2))))});
    Function g_disc_H3("g_disc_H", {xx3, uu3, dt3, lambda_disc3, p},
                        {hessian(dot(lambda_disc3, g_disc3(arg_g_disc3)[0]),
                                    vertcat(vertcat(xx3(Slice(),0),
                                                    uu3(Slice(),0),
                                                    xx3(Slice(),1)),
                                    vertcat(uu3(Slice(),1),
                                            xx3(Slice(),2))))});
    myBlocks.add_g_disc({g_disc2, g_disc3}, {g_disc_J2, g_disc_J3}, 
                        {g_disc_H2, g_disc_H3});

    // Initialization
    SX t = SX::sym("t", 1, 1);
    double T = 5.0;
    Function x_init("x_init", {t}, {vertcat(vertcat(3*t/T, 0*t/T),
                                                       SX::zeros(2,1))});
    Function u_init("u_init", {t}, {SX::zeros(2, 1)});
    myBlocks.add_inits(x_init, u_init);
    // cout<<"inits added"<<endl;

    return myBlocks;
};

BuildingBlocks makeBlocksManyObstacles(){
	int nx = 4;
    int nu = 2;
    double alpha = 10.0;

    // Helper variables
    SX xk = SX::sym("xk", nx, 1);
    SX uk = SX::sym("uk", nu, 1);
	SX t = SX::sym("t", 1, 1);
    SX p = SX::sym("p", 0, 0);
    std::vector<SX> arg = {xk, t, p};

    BuildingBlocks myBlocks = BuildingBlocks();
    myBlocks.add_sizes(nx, nu, true);

    // Initial cost
    Function Phi_0("Phi_0", arg, {alpha*t});
    Function Phi_0_J("Phi_0_J", arg, {jacobian(Phi_0(arg)[0], vertcat(t, xk))});
    Function Phi_0_H("Phi_0_H", arg, {hessian(Phi_0(arg)[0], vertcat(t, xk))});
    myBlocks.add_Phi_0(Phi_0, Phi_0_J, Phi_0_H);

    // Final cost
    Function Phi_f("Phi_f", arg, {0.0});
    Function Phi_f_J("Phi_f_J", arg, {jacobian(Phi_f(arg)[0], xk)});
    Function Phi_f_H("Phi_f_H", arg, {hessian(Phi_f(arg)[0], xk)});
    myBlocks.add_Phi_f(Phi_f, Phi_f_J, Phi_f_H);

    // Stage cost
	SX dt = SX::sym("dt", 1, 1);
    std::vector<SX> arg_phi = {xk, uk, dt, t, p};
    Function phi("phi", arg_phi, {0.0});
    Function phi_J("phi_J", arg_phi, {jacobian(phi(arg_phi)[0], 
                                              vertcat(t, xk, uk))});
    Function phi_H("phi_H", arg_phi, {hessian(phi(arg_phi)[0], 
                                              vertcat(t, xk, uk))});
    myBlocks.add_phi(phi, phi_J, phi_H);

    // Initial constraint
    SX lambda_0 = SX::sym("lambda_0", 5, 1);
    SX p_g0 = SX::sym("p_g0", 4, 1);
    SX init_cond = SX::zeros(4, 1);
    init_cond(0) = xk(0) - p_g0(0); init_cond(1) = xk(1) - p_g0(1);
    init_cond(2) = xk(2) - p_g0(2); init_cond(3) = xk(3) - p_g0(3);
    std::vector<SX> arg_g0 = {xk, t, p_g0};

    Function g0("g0", arg_g0, {vertcat(init_cond, t)});
    std::vector<double> g0_lb = {0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> g0_ub = {0.0, 0.0, 0.0, 0.0, 1000.0};
    Function g0_J("g0_J", arg_g0, {jacobian(g0(arg_g0)[0], vertcat(t, xk))});
    Function g0_H("g0_H", {xk, t, lambda_0, p},
                  {hessian(dot(lambda_0, g0(arg_g0)[0]), vertcat(t, xk))});
    myBlocks.add_g0(g0, g0_J, g0_H, g0_lb, g0_ub);

    // Final constraint
    SX lambda_T = SX::sym("lambda_T", 3, 1);
    SX p_gT = SX::sym("p_gT", 3, 1);
    std::vector<SX> arg_gT = {xk, p_gT};
    Function gT("gT", arg_gT, {vertcat(vertcat(xk(0)-p_gT(0),xk(1)-p_gT(1)),
                                               xk(2)-p_gT(2))});
    std::vector<double> gT_lb = {0.0, 0.0, 0.0};
    std::vector<double> gT_ub = {0.0, 0.0, 0.0};
    Function gT_J("gT_J", arg_gT, {jacobian(gT(arg_gT)[0], xk)});
    Function gT_H("gT_H", {xk, lambda_T, p_gT},
                  {hessian(dot(lambda_T, gT(arg_gT)[0]), xk)});
    myBlocks.add_gT(gT, gT_J, gT_H, gT_lb, gT_ub, 3);

    // Fixed constraints
    Function g_fixed("g_fixed", {xk, uk, p}, {vertcat(xk(2), xk(3), uk)});
    std::vector<double> g_fixed_lb = {0.0, -M_PI/2.0, -1.0, -M_PI/3.0};
    std::vector<double> g_fixed_ub = {2.0, M_PI/2.0, 1.0, M_PI/3.0};
    int nb_g_fixed_ = g_fixed.size1_out(0);
    SX lambda_g_fixed = SX::sym("lambda_g_fixed", nb_g_fixed_, 1);
    std::vector<SX> arg_g_fixed = {xk, uk, p};
    Function g_fixed_J("g_fixed_J", {xk, uk, p},
                                  {jacobian(g_fixed(arg_g_fixed)[0],
                                            vertcat(xk, uk))});
    Function g_fixed_H("g_fixed_H", {xk, uk, lambda_g_fixed, p},
                                  {hessian(dot(lambda_g_fixed,
                                               g_fixed(arg_g_fixed)[0]),
                                           vertcat(xk, uk))});
    myBlocks.add_g_fixed(g_fixed, g_fixed_J, g_fixed_H, g_fixed_lb,
                         g_fixed_ub);

    //  2: obstacle-avoidance constraint
    SX p_extra_2 = SX::sym("p_extra_2", 3, 1);
    Function g_extra_2("g_extra_2", {xk, uk, p_extra_2},
                                  {pow(xk(0)-p_extra_2(0), 2) +
                                   pow(xk(1)-p_extra_2(1), 2) -
                                   pow(p_extra_2(2), 2)});
    std::vector<double> lb_2 = {0.0};
    std::vector<double> ub_2 = {1.0e10};
    SX lambda_g_extra_2 = SX::sym("lambda_extra_2", g_extra_2.size1_out(0), 1);
    std::vector<SX> arg_g_extra_2 = {xk, uk, p_extra_2};
    Function g_extra_J_2("g_extra_J_2", {xk, uk, p_extra_2},
                                    {jacobian(g_extra_2(arg_g_extra_2)[0],
                                              vertcat(xk, uk))});
    Function g_extra_H_2("g_extra_H_2", {xk, uk, lambda_g_extra_2, p_extra_2},
                                    {hessian(
                                        dot(lambda_g_extra_2,
                                            g_extra_2(arg_g_extra_2)[0]),
                                        vertcat(xk, uk))});

    myBlocks.add_g_extra({g_extra_2}, {g_extra_J_2}, {g_extra_H_2}, {lb_2}, 
						 {ub_2}, {3});

    // Discretization constraints
	SX lambda_disc = SX::sym("lambda_disc", 8, 1);
	SX dts = SX::sym("dts", 2, 1);
	SX xx = SX::sym("xx", nx, 3);
	SX uu = SX::sym("uu", nu, 2);
	SX rhs1 = vertcat(vertcat(xx(2,0)*cos(xx(3,0)),
								xx(2,0)*sin(xx(3,0))),
						vertcat(uu(0,0), xx(2,0)*tan(uu(1,0))));
	SX rhs2 = vertcat(vertcat(xx(2,1)*cos(xx(3,1)),
								xx(2,1)*sin(xx(3,1))),
						vertcat(uu(0,1), xx(2,1)*tan(uu(1,1))));
	const std::vector<SX> arg_g_disc = {xx, uu, dts, t, p};
	Function g_disc("g_disc", {xx, uu, dts, t, p},
								{vertcat(xx(Slice(),1) -
										(xx(Slice(),0) + t*dts(0)*rhs1),
										xx(Slice(),2) -
										(xx(Slice(),1) + t*dts(1)*rhs2))});
	Function g_disc_J("g_disc_J", {xx, uu, dts, t, p},
									{jacobian(g_disc(arg_g_disc)[0], 
											vertcat(t,
													vertcat(xx(Slice(),0),
															uu(Slice(),0),
															xx(Slice(),1)),
													vertcat(uu(Slice(),1),
															xx(Slice(),2))))});
	Function g_disc_H("g_disc_H", {xx, uu, dts, t, lambda_disc, p},
						{hessian(dot(lambda_disc, g_disc(arg_g_disc)[0]),
									vertcat(t,
											vertcat(xx(Slice(),0),
													uu(Slice(),0),
													xx(Slice(),1)),
											vertcat(uu(Slice(),1),
													xx(Slice(),2))))});
	myBlocks.add_g_disc({g_disc}, {g_disc_J}, {g_disc_H});

    // Initialization
	double T = 1.0;
    Function x_init("x_init", {t}, {vertcat(3*t, SX::zeros(3,1))});
    Function u_init("u_init", {t}, {SX::zeros(2,1)});
    myBlocks.add_inits(x_init, u_init, 3.0);

    return myBlocks;
};

BuildingBlocks makeBlocksMoonlander(){
    int nx = 2;
    int nu = 1;

    BuildingBlocks myBlocks = BuildingBlocks();
    myBlocks.add_sizes(nx, nu, true);

    // Helper variables
    SX xk = SX::sym("xk", nx, 1);
    SX uk = SX::sym("uk", nu, 1);
    SX t = SX::sym("t", 1, 1);
    SX p = SX::sym("p", 0, 0);
    std::vector<SX> arg = {xk, t, p};

    // Initial cost
    Function Phi_0("Phi_0", arg, {t});
    Function Phi_0_J("Phi_0_J", arg, 
                     {jacobian(Phi_0(arg)[0], vertcat(t, xk))});
    Function Phi_0_H("Phi_0_H", arg, {hessian(Phi_0(arg)[0], vertcat(t, xk))});
    myBlocks.add_Phi_0(Phi_0, Phi_0_J, Phi_0_H);

    // Final cost
    Function Phi_f("Phi_f", arg, {0.0});
    Function Phi_f_J("Phi_f_J", arg, 
                     {jacobian(Phi_f(arg)[0], vertcat(t, xk))});
    Function Phi_f_H("Phi_f_H", arg, {hessian(Phi_f(arg)[0], vertcat(t, xk))});
    myBlocks.add_Phi_f(Phi_f, Phi_f_J, Phi_f_H);

    // Stage cost
    SX dt = SX::sym("dt", 1, 1);
    std::vector<SX> arg_phi = {xk, uk, dt, t, p};
    Function phi("phi", arg_phi, {0.0});
    Function phi_J("phi_J", arg_phi,
                    {jacobian(phi(arg_phi)[0], vertcat(t, xk, uk))});
    Function phi_H("phi_H", arg_phi, 
                    {hessian(phi(arg_phi)[0], vertcat(t, xk, uk))});
    myBlocks.add_phi(phi, phi_J, phi_H);

    // Initial constraint
    SX lambda_0 = SX::sym("lambda_0", 3, 1);
    std::vector<SX> arg_g0 = {xk, t, p};
    Function g0("g0", {xk, t, p}, {vertcat(xk(0), xk(1), t)});
    std::vector<double> g0_lb = {1.0, 0.0, 0.0};
    std::vector<double> g0_ub = {1.0, 0.0, 20};
    Function g0_J("g0_J", {xk, t, p}, 
                             {jacobian(g0(arg_g0)[0], vertcat(t, xk))});
    Function g0_H("g0_H", {xk, t, lambda_0, p},
                             {hessian(dot(lambda_0, g0(arg_g0)[0]), 
                                      vertcat(t, xk))});
    myBlocks.add_g0(g0, g0_J, g0_H, g0_lb, g0_ub);

    // Final constraint
    SX lambda_T = SX::sym("lambda_T", 2, 1);
    std::vector<SX> arg_gT = {xk, p};
    Function gT("gT", {xk, p}, {vertcat(xk(0), xk(1))});
    std::vector<double> gT_lb = {0.0, 0.0};
    std::vector<double> gT_ub = {0.0, 0.0};
    Function gT_J("gT_J", {xk, p}, {jacobian(gT(arg_gT)[0], xk)});
    Function gT_H("gT_H", {xk, lambda_T, p},
                             {hessian(dot(lambda_T, gT(arg_gT)[0]), xk)});
    myBlocks.add_gT(gT, gT_J, gT_H, gT_lb, gT_ub);

    // Fixed constraints
    Function g_fixed("g_fixed", {xk, uk, p}, {vertcat(xk(0), uk)});
    std::vector<double> g_fixed_lb = {0.0, 0.0};
    std::vector<double> g_fixed_ub = {1e38, 8.0};
    int nb_g_fixed_ = g_fixed.size1_out(0);
    SX lambda_g_fixed = SX::sym("lambda_g_fixed", nb_g_fixed_, 1);
    std::vector<SX> arg_g_fixed = {xk, uk, p};
    Function g_fixed_J("g_fixed_J", {xk, uk, p},
                                  {jacobian(g_fixed(arg_g_fixed)[0],
                                            vertcat(xk, uk))});
    Function g_fixed_H("g_fixed_H", {xk, uk, lambda_g_fixed, p},
                                  {hessian(dot(lambda_g_fixed,
                                               g_fixed(arg_g_fixed)[0]),
                                           vertcat(xk, uk))});
    myBlocks.add_g_fixed(g_fixed, g_fixed_J, g_fixed_H, g_fixed_lb,
                         g_fixed_ub);

    myBlocks.add_g_extra({}, {}, {}, {}, {});

    // Discretization constraints
    int nb_steps_min = 3;
    int nb_steps_max = 8;
    std::vector<Function> g_discs(nb_steps_max - nb_steps_min + 1);
    std::vector<Function> g_disc_Js(nb_steps_max - nb_steps_min + 1);
    std::vector<Function> g_disc_Hs(nb_steps_max - nb_steps_min + 1);

    Function rhs("rhs", {t, xk, uk}, {vertcat(t*xk(1), t*(-3.8 + uk))});
    std::vector<SX> arg_rhs(3);
    arg_rhs[0] = t;

    for (int nb_steps = nb_steps_min; nb_steps <= nb_steps_max; nb_steps++){
        SX lambda_disc = SX::sym("lambda_disc", nx*(nb_steps-1), 1);
        SX dts = SX::sym("dt", nb_steps-1, 1);
        SX xx = SX::sym("xx", nx, nb_steps);
        SX uu = SX::sym("uu", nu, nb_steps-1);
        SX vars = SX::zeros(1 + xx.size1()*xx.size2() + uu.size1()*uu.size2());
        vars(0) = t;
        int vars_pointer = 1;
        for (int i = 0; i < nb_steps-1; i++){
            vars(Slice(vars_pointer, vars_pointer + nx)) = xx(Slice(), i);
            vars_pointer += nx;
            vars(Slice(vars_pointer, vars_pointer + nu)) = uu(Slice(), i);
            vars_pointer += nu;
        }
        vars(Slice(vars_pointer, vars_pointer + nx)) = xx(Slice(), nb_steps-1);

        // Construct collocation grid
        SX interval_length = SX::zeros(1,1);
        for (int i = 0; i < nb_steps-1; i++){
            interval_length += dts(i);
        }
        std::vector<double> tau = collocation_points(nb_steps-1);
        DM D; DM C; DM B;
        collocation_coeff(tau, C, D, B);
        SX Cscaled = SX(C);
        for (int i = 0; i < C.size1(); i++){
            for (int j = 0; j < C.size2(); j++){
                Cscaled(i,j) = SX(C(i,j))/(interval_length);
            }
        }

        // compute derivatives of all states at the collocation points
        arg_rhs[1] = xx(Slice(), Slice(1, xx.size2()));
        arg_rhs[2] = uu;
        SX output_col = reshape(transpose(rhs(arg_rhs)[0]), nx*(nb_steps-1), 1);

        // compute derivatives of the interpolant at all collocation points 
        // for all states
        SX ders_col = reshape(mtimes(transpose(Cscaled), transpose(xx)), 
                              nx*(nb_steps-1), 1);

        const std::vector<SX> arg_g_disc = {xx, uu, dts, t, p};
        Function g_disc("g_disc", {xx, uu, dts, t, p},
                                    {ders_col - output_col});
                                    // {ders_col});
        Function g_disc_J("g_disc_J", {xx, uu, dts, t, p},
                                     {jacobian(g_disc(arg_g_disc)[0], vars)});
        Function g_disc_H("g_disc_H", {xx, uu, dts, t, lambda_disc, p},
                            {hessian(dot(lambda_disc, g_disc(arg_g_disc)[0]),
                                     vars)});
        g_discs[nb_steps-nb_steps_min] = g_disc;
        g_disc_Js[nb_steps-nb_steps_min] = g_disc_J;
        g_disc_Hs[nb_steps-nb_steps_min] = g_disc_H;
    }
    myBlocks.add_g_disc(g_discs, g_disc_Js, g_disc_Hs);

    // Initialization
    double T = 5.0;
    Function x_init("x_init", {t}, {SX::zeros(2,1)});
    Function u_init("u_init", {t}, {0.0});
    myBlocks.add_inits(x_init, u_init, 5.0);

    return myBlocks;
};