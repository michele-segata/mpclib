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

#include "MPCProblem.h"
#ifdef ENABLE_DISCRETIZATION
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>
#endif

void MPCProblem::print_message(const std::string &msg) {
    if (debug)
        std::cerr << msg;
}

#ifdef ENABLE_DISCRETIZATION

double MPCProblem::exp_A_t(double t, void *p) {
    // value to be returned
    double v;
    // get integration parameters
    struct integrate_params *params = (struct integrate_params *)p;
    // create a copy of the A matrix (we can't modify it)
    gsl_matrix *M = gsl_matrix_alloc(params->A->size1, params->A->size2);
    // matrix where we store the matrix exponential
    gsl_matrix *exp_M = gsl_matrix_alloc(params->A->size1, params->A->size2);

    // copy values of A into M
    gsl_matrix_memcpy(M, params->A);
    // compute M * x
    gsl_matrix_scale(M, t);
    // perform matrix exponentiation
    gsl_linalg_exponential_ss(M, exp_M, 0);
    // get only the value we care about
    v = gsl_matrix_get(exp_M, params->r, params->c);

    gsl_matrix_free(exp_M);
    gsl_matrix_free(M);

    return v;
}

bool MPCProblem::discretize_state_space(const Matrix<double> &A,
                                        const Matrix<double> &B,
                                        Matrix<double> &Ad, Matrix<double> &Bd,
                                        double dt) {

    // variables for gsl integration
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (10000);
    gsl_function f;
    double iv, err;
    // convert A to a gsl matrix for computation
    gsl_matrix *gsl_A = to_gsl(A);
    gsl_matrix *gsl_scaled_A = to_gsl(A);
    // discretized A matrix
    gsl_matrix *gsl_Ad = gsl_matrix_alloc(gsl_A->size1, gsl_A->size2);
    // result of the integral over the A matrix
    gsl_matrix *gsl_Ai = gsl_matrix_alloc(gsl_A->size1, gsl_A->size2);
    Matrix<double> Ai;

    // compute A * dt
    gsl_matrix_scale(gsl_scaled_A, dt);
    // perform matrix exponentiation to obtain discretized A matrix
    gsl_linalg_exponential_ss(gsl_scaled_A, gsl_Ad, 0);
    // copy result to return matrix
    Ad = from_gsl(gsl_Ad);

    // set function to integrate
    f.function = exp_A_t;
    // set integration parameters
    struct integrate_params p;
    f.params = &p;
    p.A = gsl_A;
    for (int r = 0; r < gsl_A->size1; r++) {
        for (int c = 0; c < gsl_A->size2; c++) {
            p.r = r;
            p.c = c;
            gsl_integration_qags(&f, 0, dt, 0, 1e-7, 10000, w, &iv, &err);
            gsl_matrix_set(gsl_Ai, r, c, iv);
        }
    }
    Ai = from_gsl(gsl_Ai);
    // multiply Ai by B to get Bd
    Bd = multiply(Ai, B);

    gsl_matrix_free(gsl_Ai);
    gsl_matrix_free(gsl_Ad);
    gsl_matrix_free(gsl_scaled_A);
    gsl_matrix_free(gsl_A);
    gsl_integration_workspace_free(w);

    return true;
}

#endif

void MPCProblem::setup_variables() {
    // x0 variable
    Variable *x0 = qp.vector_variable(n, "x0");
    qp.add_variable(x0);
    this->x0 = Var(x0);
    // error variable
    Variable *e_k = qp.vector_variable(q, "e_k");
    Variable *e = qp.vector_variable(T+1, e_k, "e");
    qp.add_variable(e);
    delete e_k;
    this->e = Var(e);
    // u variable, which includes u_0 as well
    Variable *u_k = qp.vector_variable(p, "u_k");
    Variable *u = qp.vector_variable(T+1, u_k, "u");
    qp.add_variable(u);
    delete u_k;
    this->u = Var(u);
    // du variable
    if (minimize_du) {
        // du variable
        Variable *du_k = qp.vector_variable(p, "du_k");
        Variable *du = qp.vector_variable(T+1, du_k, "du");
        qp.add_variable(du);
        delete du_k;
        this->du = Var(du);
    }
    if (output_slack) {
        // slack on output variable
        Variable *eps_y = qp.vector_variable(1, "eps_y");
        qp.add_variable(eps_y);
        this->eps_y = Var(eps_y);
    }
    if (control_slack) {
        // slack on control variable
        Variable *eps_u = qp.vector_variable(1, "eps_u");
        qp.add_variable(eps_u);
        this->eps_u = Var(eps_u);
    }
    if (control_derivative_slack) {
        // slack on control derivative variable
        Variable *eps_du = qp.vector_variable(1, "eps_du");
        qp.add_variable(eps_du);
        this->eps_du = Var(eps_du);
    }
    if (terminal_constraint && terminal_slack) {
        // slack on terminal constraint
        Variable *eps_t = qp.vector_variable(1, "eps_t");
        qp.add_variable(eps_t);
        this->eps_t = Var(eps_t);
    }
    s = qp.get_main_variable();
}

void MPCProblem::setup_problem() {
    // get the needed powers of the A matrix
    std::vector<Matrix<double> > Ak = get_powers(A, T);
    // precompute matrix multiplications to save time
    std::vector<Matrix<double> > C1Ak, C1AkB, nC1Ak, nC1AkB;
    std::vector<Matrix<double> > C2Ak, C2AkB, nC2Ak, nC2AkB;
    for (int i = 0; i <= T; i++) {
        C1Ak.push_back(multiply(C1, Ak[i]));
        nC1Ak.push_back(-C1Ak[C1Ak.size()-1]);
        C1AkB.push_back(multiply(C1Ak[C1Ak.size()-1], B));
        nC1AkB.push_back(-C1AkB[C1AkB.size()-1]);
        if (y_max.size() != 0 || y_min.size() != 0) {
            C2Ak.push_back(multiply(C2, Ak[i]));
            nC2Ak.push_back(-C2Ak[C2Ak.size() - 1]);
            C2AkB.push_back(multiply(C2Ak[C2Ak.size() - 1], B));
            nC2AkB.push_back(-C2AkB[C2AkB.size() - 1]);
        }
    }
    // setup some commonly used matrices and vectors
    Matrix<double> Ip = identity(p);
    Matrix<double> nIp = -Ip;
    Matrix<double> nIpTs = nIp * ts;
    Vector<double> zp = zero_vector(p);
    Matrix<double> Iq = identity(q);
    Vector<double> zq = zero_vector(q);
    Matrix<double> y_slack_M(1, q2, 1);
    Matrix<double> n_y_slack_M = -y_slack_M;
    Matrix<double> u_slack_M(1, p, 1);
    Matrix<double> n_u_slack_M = -u_slack_M;
    Matrix<double> t_slack_M(1, q, 1);
    Matrix<double> n_t_slack_M = -t_slack_M;
    Vector<double> ref(q);
    // add initial state constraint: x0 = x*
    init_state.set_variable(s);
    init_state.set_constraint_variable(x0, identity(n));
    init_state.set_known_term(init_x);
    qp.add_equality_constraint(init_state);
    // add initial control constraint: u0 = u*
    init_control.set_variable(s);
    init_control.set_constraint_variable(u[0], Ip);
    init_control.set_known_term(init_u);
    qp.add_equality_constraint(init_control);
    // error constraints e_k = C1 * A^k * x0 + ...
    for (int k = 0; k <= T; k++) {
        Constraint error(s);
        error.set_constraint_variable(e[k], Iq);
        error.set_constraint_variable(x0, nC1Ak[k]);
        for (int i = 0; i <= k - 1; i++)
            error.set_constraint_variable(u[i], nC1AkB[k-i-1]);
        if (r.size() == q) {
            // user specified a single reference, so we simply repeat it
            error.set_known_term(-r);
        }
        else {
            // if there is a signal to track, then copy the part of the signal
            // we care about in this step
            for (int i = 0; i < q; i++)
                ref[i] = -r[k * q + i];
            error.set_known_term(ref);
        }
        qp.add_equality_constraint(error);
        // save the constraint to be able to modify the reference vector
        // afterwards
        this->error.push_back(error);
    }
    // constraint u[k+1] = u[k] + du[k] * ts
    if (minimize_du) {
        for (int k = 0; k <= T - 1; k++) {
            Constraint control(s);
            control.set_constraint_variable(u[k + 1], Ip);
            control.set_constraint_variable(u[k], nIp);
            control.set_constraint_variable(du[k], nIpTs);
            control.set_known_term(zp);
            qp.add_equality_constraint(control);
        }
    }
    // terminal condition constraint: e_T = 0
    if (terminal_constraint) {
        if (!terminal_slack) {
            Constraint terminal(s);
            terminal.set_constraint_variable(e[T], Iq);
            terminal.set_known_term(zero_vector(q));
            qp.add_equality_constraint(terminal);
        }
        else {
            // -eps_t <= e_T <= eps_t
            // e_T - eps_t <= 0
            Constraint terminal1(s);
            terminal1.set_constraint_variable(e[T], Iq);
            terminal1.set_constraint_variable(eps_t, n_t_slack_M);
            terminal1.set_known_term(zero_vector(q));
            qp.add_leq_constraint(terminal1);
            // e_T + eps_t >= 0
            Constraint terminal2(s);
            terminal2.set_constraint_variable(e[T], Iq);
            terminal2.set_constraint_variable(eps_t, t_slack_M);
            terminal2.set_known_term(zero_vector(q));
            qp.add_geq_constraint(terminal2);
        }
    }

    // control constraints
    // u_k >= u_min - eps_u
    if (u_min.size() != 0) {
        for (int k = 0; k <= T - 1; k++) {
            Constraint u_min_c(s);
            u_min_c.set_constraint_variable(u[k + 1], Ip);
            if (control_slack)
                u_min_c.set_constraint_variable(eps_u, u_slack_M);
            u_min_c.set_known_term(u_min);
            qp.add_geq_constraint(u_min_c);
            c_u_min.push_back(u_min_c);
        }
    }
    // u_k <= u_max + eps_u
    if (u_max.size() != 0) {
        for (int k = 0; k <= T - 1; k++) {
            Constraint u_max_c(s);
            u_max_c.set_constraint_variable(u[k + 1], Ip);
            if (control_slack)
                u_max_c.set_constraint_variable(eps_u, n_u_slack_M);
            u_max_c.set_known_term(u_max);
            qp.add_leq_constraint(u_max_c);
            c_u_max.push_back(u_max_c);
        }
    }
    // u_k+1 - u_k >= du_min + eps_du
    if (du_min.size() != 0) {
        // if we don't minimize by du but we have du_min constraint, we add
        // u_0 - u* >= du_min => u_0 >= du_min + u*
        if (!minimize_du) {
            Constraint du_min_c(s);
            du_min_c.set_constraint_variable(u[0], Ip);
            if (control_derivative_slack)
                du_min_c.set_constraint_variable(eps_du, u_slack_M);
            du_min_c.set_known_term(du_min + init_u);
            qp.add_geq_constraint(du_min_c);
            c_du_min.push_back(du_min_c);
        }
        for (int k = 0; k <= T - 1; k++) {
            Constraint du_min_c(s);
            du_min_c.set_constraint_variable(u[k + 1], Ip);
            du_min_c.set_constraint_variable(u[k], nIp);
            if (control_derivative_slack)
                du_min_c.set_constraint_variable(eps_du, u_slack_M);
            du_min_c.set_known_term(du_min);
            qp.add_geq_constraint(du_min_c);
            c_du_min.push_back(du_min_c);
        }
    }
    // u_k+1 - u_k <= du_max + eps_du
    if (du_max.size() != 0) {
        // if we don't minimize by du but we have du_max constraint, we add
        // u_0 - u* <= du_max => u_0 <= du_max + u*
        if (!minimize_du) {
            Constraint du_max_c(s);
            du_max_c.set_constraint_variable(u[0], Ip);
            if (control_derivative_slack)
                du_max_c.set_constraint_variable(eps_du, n_u_slack_M);
            du_max_c.set_known_term(du_max + init_u);
            qp.add_leq_constraint(du_max_c);
            c_du_max.push_back(du_max_c);
        }
        for (int k = 0; k < T - 1; k++) {
            Constraint du_max_c(s);
            du_max_c.set_constraint_variable(u[k + 1], Ip);
            du_max_c.set_constraint_variable(u[k], nIp);
            if (control_derivative_slack)
                du_max_c.set_constraint_variable(eps_du, n_u_slack_M);
            du_max_c.set_known_term(du_max);
            qp.add_leq_constraint(du_max_c);
            c_du_max.push_back(du_max_c);
        }
    }
    // output constraints
    // C2 * x_k <= y_max, i.e., C2 * A^k+1 * x0 + ... <= y_max + eps_y
    if (y_max.size() != 0) {
        for (int k = 0; k <= T - 1; k++) {
            Constraint y_max_c(s);
            y_max_c.set_constraint_variable(x0, C2Ak[k+1]);
            for (int i = 0; i <= k; i++) {
                y_max_c.set_constraint_variable(u[i], C2AkB[k-i]);
                if (output_slack)
                    y_max_c.set_constraint_variable(eps_y, n_y_slack_M);
            }
            y_max_c.set_known_term(y_max);
            qp.add_leq_constraint(y_max_c);
            c_y_max.push_back(y_max_c);
        }
    }
    // C2 * x_k >= y_min, i.e., C2 * A^k+1 * x0 + ... >= y_min - eps_y
    if (y_min.size() != 0) {
        for (int k = 0; k <= T - 1; k++) {
            Constraint y_min_c(s);
            y_min_c.set_constraint_variable(x0, C2Ak[k+1]);
            for (int i = 0; i <= k; i++) {
                y_min_c.set_constraint_variable(u[i], C2AkB[k-i]);
                if (output_slack)
                    y_min_c.set_constraint_variable(eps_y, y_slack_M);
            }
            y_min_c.set_known_term(y_min);
            qp.add_geq_constraint(y_min_c);
            c_y_min.push_back(y_min_c);
        }
    }
    // if used, slack variables must be greater than 0
    if (control_slack) {
        Constraint u_slack(s);
        u_slack.set_constraint_variable(eps_u, 1);
        u_slack.set_known_term(0);
        qp.add_geq_constraint(u_slack);
    }
    if (control_derivative_slack) {
        Constraint du_slack(s);
        du_slack.set_constraint_variable(eps_du, 1);
        du_slack.set_known_term(0);
        qp.add_geq_constraint(du_slack);
    }
    if (output_slack) {
        Constraint out_slack(s);
        out_slack.set_constraint_variable(eps_y, 1);
        out_slack.set_known_term(0);
        qp.add_geq_constraint(out_slack);
    }
}

bool MPCProblem::set_state_space_matrices(const Matrix<double> &A,
                                          const Matrix<double> &B,
                                          const Matrix<double> &C1,
                                          const Matrix<double> &C2) {
    if (A.nrows() != n || A.ncols() != n) {
        print_message("Matrix A dimensions are not compatible with the size "
                      "of the state variable\n");
        return false;
    }
    if (B.nrows() != n || B.ncols() != p) {
        print_message("Matrix B dimensions are not compatible with the size "
                      "of the state variable and of the control vector\n");
        return false;
    }
    if (C1.nrows() != q || C1.ncols() != n) {
        print_message("Matrix C1 dimensions are not compatible with the size "
                      "of the state variable and of the output vector\n");
        return false;
    }
    if (C2.nrows() != q2 || C2.ncols() != n) {
        print_message("Matrix C2 dimensions are not compatible with the size "
                      "of the state variable and of the output constraint "
                      "vector\n");
        return false;
    }
    this->A = A;
    this->B = B;
    this->C1 = C1;
    this->C2 = C2;
    return true;
}

bool MPCProblem::set_initial_state(const Vector<double> &x0) {
    if (x0.size() != n) {
        print_message("x0 length is not compatible with the size of the "
                      "state vector\n");
        return false;
    }
    init_x = x0;
    return true;
}

bool MPCProblem::update_initial_state(const Vector<double> &x0) {
    if (x0.size() != n) {
        print_message("x0 length is not compatible with the size of the "
                      "state vector\n");
        return false;
    }
    init_x = x0;
    init_state.set_known_term(init_x);
    qp.update_equality_constraint(init_state, false, true);
    return true;
}

bool MPCProblem::set_initial_control(const Vector<double> &u0) {
    if (u0.size() != p) {
        print_message("u0 length is not compatible with the size of the "
                      "control vector\n");
        return false;
    }
    init_u = u0;
    return true;
}

bool MPCProblem::update_initial_control(const Vector<double> &u0) {
    if (u0.size() != p) {
        print_message("u0 length is not compatible with the size of the "
                      "state vector\n");
        return false;
    }
    init_u = u0;
    init_control.set_known_term(init_u);
    qp.update_equality_constraint(init_control, false, true);
    return true;
}

bool MPCProblem::set_reference_vector(const Vector<double> &r) {
    if (r.size() != q && r.size() < q * (T+1)) {
        print_message("r length is not compatible with the size of the "
                      "output vector\n");
        return false;
    }
    this->r = r;
    return true;
}

bool MPCProblem::update_reference_vector(const Vector<double> &r, int index) {
    if (r.size() != q) {
        print_message("r length is not compatible with the size of the "
                      "output vector\n");
        return false;
    }
    if (index == -1) {
        for (int k = 0; k <= T; k++) {
            error[k].set_known_term(-r);
            qp.update_equality_constraint(error[k], false, true);
        }
    }
    else {
        if (index >= error.size())
            return false;
        error[index].set_known_term(-r);
        qp.update_equality_constraint(error[index], false, true);
    }
    return true;
}

bool MPCProblem::seq_Q_matrix(const Matrix<double> &Q) {
    if (Q.nrows() != s.get_size() || Q.ncols() != s.get_size()) {
        print_message("Matrix Q dimensions are not compatible with the size "
                      "of the state variable\n");
        return false;
    }
    qp.set_Q_matrix(Q);
    return true;
}

bool MPCProblem::set_u_max(const Vector<double> &u_max) {
    if (u_max.size() != p) {
        print_message("u_max length is not compatible with the size of the "
                      "control vector\n");
        return false;
    }
    this->u_max = u_max;
    return true;
}

bool MPCProblem::set_u_min(const Vector<double> &u_min) {
    if (u_min.size() != p) {
        print_message("u_min length is not compatible with the size of the "
                      "control vector\n");
        return false;
    }
    this->u_min = u_min;
    return true;
}

bool MPCProblem::set_du_max(const Vector<double> &du_max) {
    if (du_max.size() != p) {
        print_message("du_max length is not compatible with the size of the "
                      "control vector\n");
        return false;
    }
    this->du_max = du_max * ts;
    return true;
}

bool MPCProblem::set_du_min(const Vector<double> &du_min) {
    if (du_min.size() != p) {
        print_message("du_min length is not compatible with the size of the "
                      "control vector\n");
        return false;
    }
    this->du_min = du_min * ts;
    return true;
}

bool MPCProblem::set_y_max(const Vector<double> &y_max) {
    if (y_max.size() != q2) {
        print_message("y_max length is not compatible with the size of the "
                      "output constraint vector\n");
        return false;
    }
    this->y_max = y_max;
    return true;
}

bool MPCProblem::set_y_min(const Vector<double> &y_min) {
    if (y_min.size() != q2) {
        print_message("y_min length is not compatible with the size of the "
                      "output constraint vector\n");
        return false;
    }
    this->y_min = y_min;
    return true;
}

bool MPCProblem::update_u_max_vector(const Vector<double> &u_max, int index) {
    if (u_max.size() != p) {
        print_message("u_max length is not compatible with the size of the "
                      "control vector\n");
        return false;
    }
    if (index == -1) {
        for (int k = 0; k <= T-1; k++) {
            c_u_max[k].set_known_term(u_max);
            qp.update_leq_constraint(c_u_max[k], false, true);
        }
    }
    else {
        if (index >= c_u_max.size())
            return false;
        c_u_max[index].set_known_term(u_max);
        qp.update_leq_constraint(c_u_max[index], false, true);
    }
    return true;
}

bool MPCProblem::update_u_min_vector(const Vector<double> &u_min, int index) {
    if (u_min.size() != p) {
        print_message("u_min length is not compatible with the size of the "
                      "control vector\n");
        return false;
    }
    if (index == -1) {
        for (int k = 0; k <= T-1; k++) {
            c_u_min[k].set_known_term(u_min);
            qp.update_geq_constraint(c_u_min[k], false, true);
        }
    }
    else {
        if (index >= c_u_min.size())
            return false;
        c_u_min[index].set_known_term(u_min);
        qp.update_geq_constraint(c_u_min[index], false, true);
    }
    return true;
}

bool MPCProblem::update_du_max_vector(const Vector<double> &du_max, int index) {
    if (du_max.size() != p) {
        print_message("du_max length is not compatible with the size of the "
                      "control vector\n");
        return false;
    }
    if (index == -1) {
        for (int k = 0; k <= T-1; k++) {
            c_du_max[k].set_known_term(du_max);
            qp.update_leq_constraint(c_du_max[k], false, true);
        }
    }
    else {
        if (index >= c_du_max.size())
            return false;
        c_du_max[index].set_known_term(du_max);
        qp.update_leq_constraint(c_du_max[index], false, true);
    }
    return true;
}

bool MPCProblem::update_du_min_vector(const Vector<double> &du_min, int index) {
    if (du_min.size() != p) {
        print_message("du_min length is not compatible with the size of the "
                      "control vector\n");
        return false;
    }
    if (index == -1) {
        for (int k = 0; k <= T-1; k++) {
            c_du_min[k].set_known_term(du_min);
            qp.update_geq_constraint(c_du_min[k], false, true);
        }
    }
    else {
        if (index >= c_du_min.size())
            return false;
        c_du_min[index].set_known_term(du_min);
        qp.update_geq_constraint(c_du_min[index], false, true);
    }
    return true;
}

bool MPCProblem::update_y_min_vector(const Vector<double> &y_min, int index) {
    if (y_min.size() != q2) {
        print_message("y_min length is not compatible with the size of the "
                      "output constraint vector\n");
        return false;
    }
    if (index == -1) {
        for (int k = 0; k <= T-1; k++) {
            c_y_min[k].set_known_term(y_min);
            qp.update_geq_constraint(c_y_min[k], false, true);
        }
    }
    else {
        if (index >= c_y_min.size())
            return false;
        c_y_min[index].set_known_term(y_min);
        qp.update_geq_constraint(c_y_min[index], false, true);
    }
    return true;
}

bool MPCProblem::update_y_max_vector(const Vector<double> &y_max, int index) {
    if (y_max.size() != q2) {
        print_message("y_max length is not compatible with the size of the "
                      "output constraint vector\n");
        return false;
    }
    if (index == -1) {
        for (int k = 0; k <= T-1; k++) {
            c_y_max[k].set_known_term(y_max);
            qp.update_leq_constraint(c_y_max[k], false, true);
        }
    }
    else {
        if (index >= c_y_max.size())
            return false;
        c_y_max[index].set_known_term(y_max);
        qp.update_leq_constraint(c_y_max[index], false, true);
    }
    return true;
}

double MPCProblem::solve_mpc(Vector<double> &arg) {
    // TODO: check whether things are set
    if (!problem_setup) {
        setup_problem();
        problem_setup = true;
        if (debug)
            qp.write_information();
    }
    return qp.solve_problem(arg);
}

bool MPCProblem::is_feasible(double result) {
    if (std::numeric_limits<double>::has_infinity)
        return !std::isinf(result);
    else
        return result < 1.0E300;
}

Matrix<double> MPCProblem::get_state_evolution(const Matrix<double> &C,
                                               const Vector<double> &solution) {
    int samples = T + 1;
    int states = C.nrows();
    Matrix<double> output(samples, states);
    Matrix<double> x(samples, n);
    Matrix<double> yk, xk;
    Vector<double> uk;

    copy_vector(x[0], init_x);
    yk = multiply(C, x[0]);
    copy_vector(output[0], yk);
    for (int k = 1; k <= T; k++) {
        // x_k = A x_k-1 + B u_k-1
        uk = subvector(solution, u[k-1].get_position(), u[k-1].get_size());
        xk = multiply(A, x[k-1]) + multiply(B, uk);
        copy_vector(x[k], xk);
        // y_k = C x_k
        yk = multiply(C, x[k]);
        copy_vector(output[k], yk);
    }
    return output;
}
