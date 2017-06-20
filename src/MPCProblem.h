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

#ifndef MPCPROBLEM_H
#define MPCPROBLEM_H

#include "Array.hh"
#include "QuadraticProblem.h"

/**
 * High level interface for generating and solving model predictive control
 * problems
 */
class MPCProblem {

private:

    // prediction horizon in number of steps
    int T;
    // size of the state matrix
    int n;
    // size of the input vector
    int p;
    // size of the output vector
    int q;
    // size of the output constraint vector
    int q2;
    // state matrix. this must be the state matrix of the discretized system
    Matrix<double> A;
    // control matrix
    Matrix<double> B;
    // output matrix
    Matrix<double> C1;
    // output constraint matrix
    Matrix<double> C2;
    // reference vector
    Vector<double> r;
    // initial state
    Vector<double> init_x;
    // initial control
    Vector<double> init_u;
    // quadratic problem solver
    QuadraticProblem qp;
    // sampling time
    double ts;
    // whether we care about minimization of the control derivative
    bool minimize_du;
    // whether we enable terminal cost or not
    bool terminal_constraint;
    // enable/disable slack variable on output constraint
    bool output_slack;
    // enable/disable slack variable on control constraint
    bool control_slack;
    // enable/disable slack variable on control derivative constraint
    bool control_derivative_slack;
    // enable/disable slack variable on terminal constraint
    bool terminal_slack;
    // indicates whether the problem setup is already done
    bool problem_setup;

    // main variable
    Var s;
    // initial state
    Var x0;
    // error variable
    Var e;
    // control variable
    Var u;
    // control derivative variable
    Var du;
    // slack variable on output
    Var eps_y;
    // slack variable on control
    Var eps_u;
    // slack variable on control derivative
    Var eps_du;
    // slack variable on terminal constraint
    Var eps_t;

    // constraint on initial state
    Constraint init_state;
    // constraint on initial control
    Constraint init_control;
    // constraints on error
    std::vector<Constraint> error;
    // limit constraints
    std::vector<Constraint> c_u_min, c_u_max, c_du_min, c_du_max, c_y_min,
                            c_y_max;

    // control, control derivative, and output limits
    Vector<double> u_min, u_max, du_min, du_max, y_min, y_max;

    // debug output
    bool debug;

    /**
     * Sets up the variables
     */
    void setup_variables();

    /**
     * If the debug flag is enabled, prints a message to stderr
     * @param msg the message to be printed
     */
    void print_message(const std::string &msg);

#ifdef ENABLE_DISCRETIZATION

    /**
     * Struct used to pass parameters to the function being integrated for
     * discretization
     */
    struct integrate_params {
        // continuous time state matrix
        const gsl_matrix *A;
        // row on which we are doing the integration
        int r;
        // column on which we are doing the integration
        int c;
    };

    /**
     * Function on which we compute the integral for discretization. To
     * discretize the control matrix we need to compute
     *
     * int_0^{dt} e^{A t} dt * B
     *
     * This function simply computes e^{A t} given A and t.
     * The result of the integral over the matrix is a matrix, but the
     * function being integrated must return a single value. To cope with this,
     * we must repeat the integration operation for all rows and columns, so
     * this function also requires the row and the column of the element that
     * the function should return
     *
     * @param t integration variable
     * @param p a pointer to a struct integrate_params, which includes the
     * continuous time state matrix A and the row and column indexes
     * @return the value of e^{A t}[r][c]
     */
    static double exp_A_t(double t, void *p);

#endif

public:

    /**
     * Constructor. Sets the main variables of the control problem
     * @param T prediction horizon in steps
     * @param n size of the state matrix and of the state vector
     * @param p size of the control vector
     * @param q2 size of the ouput constraint vector. if no output constraints
     * are given, this value is simply ignored
     * @param q size of the ouput vector
     * @param ts sampling time
     * @param minimize_du if set to true, the control derivative is taken
     * into account in the minimization
     * @param terminal_constraint enable/disable the terminal constraint,
     * that is e[T-1] = 0
     * @param output_slack enable/disable slack variable on output constraints
     * @param control_slack enable/disable slack variable on control
     * constraints
     * @param control_derivative_slack enable/disable slack variable on
     * control derivative constraints
     * @param terminal_slack enable/disable slack variable on the terminal
     * constraint. If the terminal constraint is disabled, this parameter is
     * ignored
     * @param debug enable/disable debug output on stderr
     */
    MPCProblem(int T, int n, int p, int q, int q2, double ts,
               bool minimize_du, bool terminal_constraint, bool output_slack,
               bool control_slack, bool control_derivative_slack,
               bool terminal_slack, bool debug = false) : s(0), x0(0), e(0),
               u(0), du(0), eps_y(0), eps_u(0), eps_du(0), eps_t(0) {
        this->T = T;
        this->n = n;
        this->p = p;
        this->q = q;
        this->q2 = q2;
        this->ts = ts;
        this->minimize_du = minimize_du;
        this->terminal_constraint = terminal_constraint;
        this->output_slack = output_slack;
        this->control_slack = control_slack;
        this->control_derivative_slack = control_derivative_slack;
        this->terminal_slack = terminal_slack;
        this->debug = debug;
        this->problem_setup = false;
        setup_variables();
    }

#ifdef ENABLE_DISCRETIZATION

    /**
     * Discretizes a continuous time state space system
     * @param A continuous time state matrix
     * @param B continuous time control matrix
     * @param Ad returned discrete time state matrix
     * @param Bd returned discrete time control matrix
     * @param dt sampling time in seconds
     * @return true if the operation succeeds, false otherwise
     */
    static bool discretize_state_space(const Matrix<double> &A,
                                       const Matrix<double> &B,
                                       Matrix<double> &Ad, Matrix<double> &Bd,
                                       double dt);

#endif


    /**
     * Sets up the quadratic programming problem to solve. This method can be
     * invoked after setting state space matrices and constraint to construct
     * problem matrices and vectors. If not called prior invoking solve_mpc(),
     * the library automatically calls the method before solving the problem.
     * Be aware that, before setup_problem() is called, it is not possible to
     * invoke update methods such as update_reference_vector().
     */
    void setup_problem();

    /**
     * Sets the A, B, and C matrices
     * @param A state matrix (discrete time)
     * @param B control matrix (discrete time)
     * @param C1 output matrix
     * @param C2 output constraint matrix. if no output constraints are given,
     * this matrix is simply ignored
     * @return false if the matrix sizes are incompatible, true otherwise
     */
    bool set_state_space_matrices(const Matrix<double> &A,
                                  const Matrix<double> &B,
                                  const Matrix<double> &C1,
                                  const Matrix<double> &C2);

    /**
     * Sets the initial state
     * @param x0 the initial state
     * @return false if the size of the given vector is incompatible, true
     * otherwise
     */
    bool set_initial_state(const Vector<double> &x0);

    /**
     * Updates the initial state, changing the value already set. This is
     * useful when solving again the same problem under different initial
     * conditions
     * @param x0 the new initial state
     * @return false if the size of the given vector is incompatible, true
     * otherwise
     */
    bool update_initial_state(const Vector<double> &x0);

    /**
     * Sets the initial control
     * @param u0 the initial control
     * @return false if the size of the given vector is incompatible, true
     * otherwise
     */
    bool set_initial_control(const Vector<double> &u0);

    /**
     * Updates the initial control, changing the value already set. This is
     * useful when solving again the same problem under different initial
     * conditions
     * @param u0 the new initial control
     * @return false if the size of the given vector is incompatible, true
     * otherwise
     */
    bool update_initial_control(const Vector<double> &u0);

    /**
     * Sets the reference value
     * @param r the reference value vector
     * @return false if the size of the given vector is incompatible, true
     * otherwise
     */
    bool set_reference_vector(const Vector<double> &r);

    /**
     * Updates the reference value
     * @param r the reference value vector
     * @param index which reference vector to modify. If set to -1, all
     * reference vectors will be modified
     * @return false if the size of the given vector is incompatible, true
     * otherwise
     */
    bool update_reference_vector(const Vector<double> &r, int index=-1);

    /**
     * Sets the Q matrix of the cost function
     * @param Q the matrix
     * @return false if the size of the given matrix is incompatible, true
     * otherwise
     */
    bool seq_Q_matrix(const Matrix<double> &Q);

    /**
     * Sets the u_max limit
     * @param u_max control limit
     * @return false if size of the given vector is incompatible, true
     * otherwise
     */
    bool set_u_max(const Vector<double> &u_max);

    /**
     * Updates the u_max constraint
     * @param u_max the constraint value vector
     * @param index which constraint vector to modify. If set to -1, all
     * constraint vectors will be modified
     * @return false if the size of the given vector is incompatible, true
     * otherwise
     */
    bool update_u_max_vector(const Vector<double> &u_max, int index=-1);

    /**
     * Sets the u_min limit
     * @param u_min control limit
     * @return false if size of the given vector is incompatible, true
     * otherwise
     */
    bool set_u_min(const Vector<double> &u_min);

    /**
     * Updates the u_min constraint
     * @param u_min the constraint value vector
     * @param index which constraint vector to modify. If set to -1, all
     * constraint vectors will be modified
     * @return false if the size of the given vector is incompatible, true
     * otherwise
     */
    bool update_u_min_vector(const Vector<double> &u_min, int index=-1);

    /**
     * Sets the du_max limit
     * @param du_max limit on control derivative. this value refers to the
     * derivative, and not the discrete step. Internally, du_max is
     * multiplied by the sampling time
     * @return false if size of the given vector is incompatible, true
     * otherwise
     */
    bool set_du_max(const Vector<double> &du_max);

    /**
     * Updates the du_max constraint
     * @param du_max the constraint value vector
     * @param index which constraint vector to modify. If set to -1, all
     * constraint vectors will be modified
     * @return false if the size of the given vector is incompatible, true
     * otherwise
     */
    bool update_du_max_vector(const Vector<double> &du_max, int index=-1);

    /**
     * Sets the du_min limit
     * @param du_min limit on control derivative. this value refers to the
     * derivative, and not the discrete step. Internally, du_min is
     * multiplied by the sampling time
     * @return false if size of the given vector is incompatible, true
     * otherwise
     */
    bool set_du_min(const Vector<double> &du_min);

    /**
     * Updates the du_min constraint
     * @param du_min the constraint value vector
     * @param index which constraint vector to modify. If set to -1, all
     * constraint vectors will be modified
     * @return false if the size of the given vector is incompatible, true
     * otherwise
     */
    bool update_du_min_vector(const Vector<double> &du_min, int index=-1);

    /**
     * Sets the y_min limit. This is interpreted as C*x_k >= y_min
     * @param y_min limit on output
     * @return false if size of the given vector is incompatible, true
     * otherwise
     */
    bool set_y_min(const Vector<double> &y_min);

    /**
     * Updates the y_min constraint
     * @param y_min the constraint value vector
     * @param index which constraint vector to modify. If set to -1, all
     * constraint vectors will be modified
     * @return false if the size of the given vector is incompatible, true
     * otherwise
     */
    bool update_y_min_vector(const Vector<double> &y_min, int index=-1);

    /**
     * Sets the y_max limit. This is interpreted as C*x_k <= y_max
     * @param y_max limit on output
     * @return false if size of the given vector is incompatible, true
     * otherwise
     */
    bool set_y_max(const Vector<double> &y_max);

    /**
     * Updates the y_max constraint
     * @param y_max the constraint value vector
     * @param index which constraint vector to modify. If set to -1, all
     * constraint vectors will be modified
     * @return false if the size of the given vector is incompatible, true
     * otherwise
     */
    bool update_y_max_vector(const Vector<double> &y_max, int index=-1);

    /**
     * Solves the mpc problem
     * @param arg the argument minimizing the cost function
     * @return the minimum value of the cost function
     */
    double solve_mpc(Vector<double> &arg);

    /**
     * Returns whether the return value of solve_mpc() is valid or if the
     * problem is infeasible.
     * @param result the minimum value computed by solve_mpc
     * @return true if the problem was feasible and the result is the minimum,
     * false if the problem was not feasible
     */
    bool is_feasible(double result);

    /**
     * Given a certain output matrix and the solution to the problem, computes
     * the evolution of the system. The function simply computes
     * y_k = C * (A * x_k + B * u_k), for k = 0..T
     * for a given output matrix C and the solution u_k computed by the solver.
     * @param C the output matrix
     * @param solution the solution computed by solve_mpc()
     * @return a matrix M[k][i] where k (the row) indicates the sample index
     * (discrete time) and i (the column) is the i-th variable of the system
     * output
     */
    Matrix<double> get_state_evolution(const Matrix<double> &C,
                                       const Vector<double> &solution);
    Var get_s() {
        return s;
    }

    Var get_x0() {
        return x0;
    }

    Var get_e() {
        return e;
    }

    Var get_u() {
        return u;
    }

    Var get_du() {
        return du;
    }

};


#endif //MPCPROBLEM_H
