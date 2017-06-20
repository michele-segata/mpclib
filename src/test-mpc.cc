//
// Copyright (c) 2017 Michele Segata <msegata@disi.unitn.it>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see http://www.gnu.org/licenses/.
//
// Tester program for an MPC problem. The problem is defined as a tracking
// problem on a vehicle where the reference speed is set to 1 and the initial
// speed and acceleration are set to 0. The plant model includes a first
// order lag on the acceleration, i.e., following the transfer function
//
// a(s) = 1 / (tau * s + 1) * u(s)
//
// where a(s) is the acceleration, u(s) is the control input, and tau is the
// lag time constant, which is set to 0.5 s. The state space representation
// of the system in continuous time is
//
// dx = A*x + B*u
// y = C1*x
//
// where A = | -1/tau, 0 |    B = | 1/tau, 0 |   C1 = | 0, 1 |
//           |      1, 0 |
//
// and x = y = | a |
//             | v |
//
// with a and v being acceleration and speed respectively. The program
// considers the discretized state space represetantion for a sampling time
// of 0.1 seconds, i.e.,
//
// x(k+1) = A * x(k) + B * u(k)
// y(k) = C1 * x(k)
//
// where k is the sample index. By discretizing A, B, and C1, we obtain
//
// A = | 0.8187 , 0 |    B = | 0.1813, 0.009365 |   C1 = | 0, 1 |
//     | 0.09063, 1 |
//
// The program writes the result on stdout in csv format. The output includes
// four fields which are the acceleration error (ea), the speed error (ev),
// the control (u) and the control derivative (du). The result can be plotted
// using the given R script. Simply run
//
// build/src/test-mpc > res.csv
// Rscript plot-mpc.R
//
// The plot is written inside mpc.pdf.

#include "Array.hh"
#include "matrix_utils.h"
#include "MPCProblem.h"
#include "tclap/CmdLine.h"
#include "tclap/ValueArg.h"
#include "tclap/SwitchArg.h"
#include "tclap/ArgException.h"
#ifdef ENABLE_CONFIG_FILE
#include "ConfigLoader.h"
#endif

using namespace std;

#ifndef ISZERO
#define ISZERO(x) (abs(x) < 1e-15)
#endif

typedef TCLAP::SwitchArg SA;
typedef TCLAP::ValueArg<double> VAD;
typedef TCLAP::ValueArg<int> VAI;
typedef TCLAP::ValueArg<string> VAS;

/**
 * Get the discrete time version of the state matrix
 * @param ts sampling time
 * @param tau actuation lag
 * @return discrete state matrix
 */
Matrix<double> get_A(double ts, double tau);

/**
 * Get the discrete time version of the control matrix
 * @param ts sampling time
 * @param tau actuation lag
 * @return discrete control matrix
 */
Matrix<double> get_B(double ts, double tau);

int main(int argc, char** argv) {

    // define command line arguments
    VAD a_min("a", "a-min", "acceleration lower bound", false, -1e9, "double");
    VAD a_max("b", "a-max", "acceleration upper bound", false, 1e9, "double");
#ifdef ENABLE_CONFIG_FILE
    VAS config("c", "config-file", "config file defining continuous time "
               "state space matrices", false, "", "string");
#endif
    VAD du_weight("d", "control-derivative-weight", "weight for control "
                  "derivative variables", false, 1.0, "double");
    VAD error_weight("e", "error-weight", "weight for speed error variables",
                     false, 1.0, "double");
    VAD u_slack("q", "u-slack", "slack variable weight for control bounds. "
                "Set to 0 for disabling it", false, 0, "double");
    VAD du_slack("r", "du-slack", "slack variable weight for control "
                 "derivative bounds. Set to 0 for disabling it", false, 0,
                 "double");
    VAD out_slack("s", "out-slack", "slack variable weight for output bounds. "
                  "Set to 0 for disabling it", false, 0, "double");
    VAI simulate("S", "simulate", "simulates the system for the specified "
                 "number of steps", false, 30, "int");
    VAD terminal_slack("f", "terminal-slack", "slack variable weight for "
                       "terminal constraint. Set to 0 for disabling it", false,
                       0, "double");
    VAD tau("t", "tau", "actuation lag in seconds", false, 0.5, "double");
    VAI horizon("T", "horizon", "time horizon in steps", false, 60, "int");
    VAD u_weight("u", "control-weight", "weight for control variables", false,
                 1.0, "double");
    VAD u_min("v", "u-min", "control lower bound", false, -1e9, "double");
    VAD u_max("w", "u-max", "control upper bound", false, 1e9, "double");
    VAD du_min("x", "du-min", "control derivative lower bound", false, -1e9,
               "double");
    VAD du_max("y", "du-max", "control derivative upper bound", false, 1e9,
               "double");
    SA no_terminal("z", "no-terminal", "disables terminal constraint. this "
                   "can be used to make the problem feasible");
    SA debug("g", "debug", "enable debug output printing information about "
             "the optimization problem");

    try {

        TCLAP::CmdLine cmd("Test application for MPCProblem class", ' ', "1.0");

        cmd.add(error_weight);
        cmd.add(u_weight);
        cmd.add(du_weight);
        cmd.add(u_min);
        cmd.add(u_max);
        cmd.add(du_min);
        cmd.add(du_max);
        cmd.add(a_min);
        cmd.add(a_max);
        cmd.add(tau);
        cmd.add(horizon);
        cmd.add(no_terminal);
        cmd.add(debug);
        cmd.add(u_slack);
        cmd.add(du_slack);
        cmd.add(out_slack);
        cmd.add(terminal_slack);
        cmd.add(simulate);
#ifdef ENABLE_CONFIG_FILE
        cmd.add(config);
#endif

        // Parse the argv array.
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException &e) {
        cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
        return 1;
    }

    //prediction horizon (in steps)
    int T = horizon.getValue();
    //number of variables in the state (state space)
    int n = 2;
    //input size
    int p = 1;
    //output size
    int q = 1;
    //output constraint size
    int q2 = 1;
    //sampling time
    double ts = 0.1;

    //state, input, output, and constraint output matrices
    Matrix<double> A, B, C1, C2;
    //reference vector
    Vector<double> ref;

#ifdef ENABLE_CONFIG_FILE
    if (config.isSet()) {
        ConfigLoader cfg(config.getValue());
        if (!cfg.parseConfig())
            return 1;
        if (!cfg.loadConfiguration())
            return 1;
        //continuous time state space matrices
        Matrix<double> Ac, Bc;
        Ac = cfg.get_A();
        Bc = cfg.get_B();
        C1 = cfg.get_C1();
        C2 = cfg.get_C2();
        ref = cfg.get_ref_vector();
        n = Ac.nrows();
        p = Bc.ncols();
        q = C1.nrows();
        q2 = C2.nrows();
        MPCProblem::discretize_state_space(Ac, Bc, A, B, ts);
    }
    else {
#endif
        A = get_A(ts, tau.getValue());
        B = get_B(ts, tau.getValue());
        C1 = get_matrix("0,1", q, n);
        C2 = get_matrix("1,0", q2, n);
        //set reference speed as 1
        ref.resize(1, q);

#ifdef ENABLE_CONFIG_FILE
    }
#endif

    //initial state
    Vector<double> init_x(0.0, n), init_u(0.0, p);

    // bounds
    Vector<double> ymin(a_min.getValue(), q2);
    Vector<double> ymax(a_max.getValue(), q2);
    Vector<double> umin(u_min.getValue(), p);
    Vector<double> umax(u_max.getValue(), p);
    Vector<double> dumin(du_min.getValue(), p);
    Vector<double> dumax(du_max.getValue(), p);

    // include control derivative in the minimization
    bool minimize_du = true;
    // enable the terminal constraint in the quadratic problem
    bool terminal_constraint = !no_terminal.getValue();
    // enable/disable slack on output limits
    bool output_slack = !ISZERO(out_slack.getValue());
    // enable/disable slack on control limits
    bool control_slack = !ISZERO(u_slack.getValue());
    // enable/disable slack on control derivative limits
    bool control_derivative_slack = !ISZERO(du_slack.getValue());
    // enable/disable slack on control derivative limits
    bool term_slack = !ISZERO(terminal_slack.getValue());

    // instantiate the problem solver
    MPCProblem mpc(T, n, p, q, q2, ts, minimize_du, terminal_constraint,
                   output_slack, control_slack, control_derivative_slack,
                   term_slack, debug.getValue());

    // set state space matrices
    mpc.set_state_space_matrices(A, B, C1, C2);
    // set control limits
    if (u_max.isSet())
        mpc.set_u_max(umax);
    if (u_min.isSet())
        mpc.set_u_min(umin);
    if (du_max.isSet())
        mpc.set_du_max(dumax);
    if (du_min.isSet())
        mpc.set_du_min(dumin);
    if (a_max.isSet())
        mpc.set_y_max(ymax);
    if (a_min.isSet())
        mpc.set_y_min(ymin);
    // set initial state, initial control, and reference
    mpc.set_initial_state(init_x);
    mpc.set_initial_control(init_u);
    mpc.set_reference_vector(ref);

    // get problem variables for setting up the weight matrix
    Var x0 = mpc.get_x0();
    Var e = mpc.get_e();
    Var u = mpc.get_u();
    Var du = mpc.get_du();

    // weight on x0. this doesn't change anything as x0 is fixed
    Vector<double> Qx(1, x0.get_size());
    // weight on speed error
    Vector<double> Qe(error_weight.getValue(), e.get_size());
    // weight on control
    Vector<double> Qu(u_weight.getValue(), u.get_size());
    // weight on control derivative
    Vector<double> Qdu;
    if (minimize_du)
        Qdu.resize(du_weight.getValue(), du.get_size());
    // weight for output slack
    Vector<double> Qout_slack(out_slack.getValue(), 1);
    // weight for control slack
    Vector<double> Qu_slack(u_slack.getValue(), 1);
    // weight for control derivative slack
    Vector<double> Qdu_slack(du_slack.getValue(), 1);
    // weight for control terminal slack
    Vector<double> Qt_slack(terminal_slack.getValue(), q);

    // merge together all weight vectors...
    Vector<double> Qv(Qx);
    Qv = bind_vectors(Qv, Qe);
    Qv = bind_vectors(Qv, Qu);
    if (minimize_du)
        Qv = bind_vectors(Qv, Qdu);
    if (output_slack)
        Qv = bind_vectors(Qv, Qout_slack);
    if (control_slack)
        Qv = bind_vectors(Qv, Qu_slack);
    if (control_derivative_slack)
        Qv = bind_vectors(Qv, Qdu_slack);
    if (term_slack)
        Qv = bind_vectors(Qv, Qt_slack);
    // and create the diagonal matrix of weights
    Matrix<double> Q = diag(Qv);
    mpc.seq_Q_matrix(Q);

    // variables including problem solution
    Vector<double> solution;
    Matrix<double> output;
    double min;

    // number of simulation steps to run. if the variable is not set, we do
    // only one
    int sim_steps = 1;
    if (simulate.isSet())
        sim_steps = simulate.getValue();

    // write the result to stdout (sim step, solution step, state, control,
    // and control derivative)

    // write csv header
    cout << "n,k,";
    for (int i = 0; i < n; i++)
        cout << "x" << i << ",";
    for (int i = 0; i < p; i++)
        cout << "u" << i << ",";
    for (int i = 0; i < p; i++) {
        cout << "du" << i;
        if (i != p - 1)
            cout << ",";
    }
    cout << "\n";

    for (int t = 0; t < sim_steps; t++) {

        // solve the problem
        min = mpc.solve_mpc(solution);
        if (!mpc.is_feasible(min)) {
            // if unfeasible, no need to continue
            break;
        }
        // get the evolution of the state given the solution
        output = mpc.get_state_evolution(identity(n), solution);

        for (int k = 0; k <= T; k++) {

            // output sim time and solution step
            cout << t << "," << k << ",";
            // output state variables
            for (int i = 0; i < output.ncols(); i++)
                cout << output[k][i] << ",";
            // output control variables
            for (int i = 0; i < u[k].get_size(); i++)
                cout << solution[u[k][i].get_position()] << ",";
            // output control derivative variables
            for (int i = 0; i < du[k].get_size(); i++) {
                if (minimize_du)
                    cout << solution[du[k][i].get_position()];
                else
                    cout << 0;

                if (i != du[k].get_size() - 1)
                    cout << ",";
            }
            cout << "\n";
        }

        // update initial state and initial control for next simulation step
        for (int i = 0; i < t; i++)
            init_x[i] = output[1][i];
        for (int i = 0; i < p; i++)
            init_u[i] = solution[u[1][i].get_position()];
        mpc.update_initial_state(init_x);
        mpc.update_initial_control(init_u);
    }

    return 0;
}

Matrix<double> get_A(double ts, double tau) {
    Matrix<double> A(2, 2);
    A[0][0] = exp(-ts/tau);
    A[0][1] = 0;
    A[1][0] = tau * (1 - exp(-ts/tau));
    A[1][1] = 1;
    return A;
}

Matrix<double> get_B(double ts, double tau) {
    Matrix<double> B(2, 1);
    B[0][0] = 1 - exp(-ts/tau);
    B[1][0] = ts + tau * (exp(-ts/tau) - 1);
    return B;
}
