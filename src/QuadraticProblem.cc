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

#include <map>
#include "QuadraticProblem.h"
#include "QuadProg++.hh"

Variable *QuadraticProblem::vector_variable(int n, const std::string &name,
                                            const std::string &base_name) {
    Variable *base = new Variable(true, base_name);
    Variable *vector = new Variable(false, name);
    vector->add_variable(base, n);
    delete base;
    return vector;
}

Variable *QuadraticProblem::vector_variable(int n,
                                            const Variable *base_variable,
                                            const std::string &name) {
    Variable *vector = new Variable(false, name);
    vector->add_variable(base_variable, n);
    return vector;
}

void QuadraticProblem::print_message(const std::string &msg) {
    if (debug)
        std::cerr << msg;
}

bool QuadraticProblem::add_variable(Variable *v) {
    VarMap::iterator i;
    if (constraint_phase) {
        print_message("Cannot add variables anymore. User already started "
                      "adding constraints\n");
        return false;
    }
    i = variables.find(v);
    if (i != variables.end()) {
        print_message("The given variable has already been added.\n");
        return false;
    }
    else {
        variables[v] = vx->get_subvariables_count();
    }
    vx->add_variable(v);
    return true;
}

Var QuadraticProblem::get_variable(const Variable *v) {
    return Var(v);
}

void QuadraticProblem::init_constraint_phase() {
    if (constraint_phase)
        return;
    constraint_phase = true;
    eq.set_variable(Var(vx));
    geq.set_variable(Var(vx));
    q0.resize(0, vx->get_size());
}

Var QuadraticProblem::get_main_variable() {
    init_constraint_phase();
    return Var(vx);
}

bool QuadraticProblem::add_equality_constraint(Constraint &c) {
    init_constraint_phase();
    if (c.get_variable().get_variable() != vx) {
        print_message("The given constraint has been built on a different "
                      "variable and cannot be added to the problem.\n");
        return false;
    }
    return eq.add_constraint(c);
}

bool QuadraticProblem::update_equality_constraint(const Constraint &c,
                                                  bool update_matrix,
                                                  bool update_vector) {
    return eq.update_constraint(c, update_matrix, update_vector);
}

bool QuadraticProblem::add_geq_constraint(Constraint &c) {
    init_constraint_phase();
    if (c.get_variable().get_variable() != vx) {
        print_message("The given constraint has been built on a different "
                      "variable and cannot be added to the problem.\n");
        return false;
    }
    return geq.add_constraint(c);
}

bool QuadraticProblem::update_geq_constraint(const Constraint &c,
                                             bool update_matrix,
                                             bool update_vector) {
    return geq.update_constraint(c, update_matrix, update_vector);
}

bool QuadraticProblem::add_leq_constraint(Constraint &c) {
    bool result;
    init_constraint_phase();
    if (c.get_variable().get_variable() != vx) {
        print_message("The given constraint has been built on a different "
                      "variable and cannot be added to the problem.\n");
        return false;
    }
    // transform <= into >=
    c.invert_constraint();
    result = geq.add_constraint(c);
    c.invert_constraint();
    return result;

}

bool QuadraticProblem::update_leq_constraint(const Constraint &c,
                                             bool update_matrix,
                                             bool update_vector) {
    // transform <= into >=
    Constraint inv(c);
    inv.invert_constraint();
    return geq.update_constraint(inv, update_matrix, update_vector);
}

bool QuadraticProblem::set_Q_matrix(const Matrix<double> &Q) {
    if (!constraint_phase) {
        print_message("Before setting the Q matrix you should enter in the "
                      "constraint phase by either adding constraints or by "
                      "obtaining the x variable via get_main_variable()\n");
        return false;
    }
    if (Q.ncols() != vx->get_size() || Q.nrows() != vx->get_size()) {
        print_message("Incompatible Q matrix. Q should be an N by N square "
                      "matrix with N being the size of the x variable\n");
        return false;
    }
    this->Q = Q;
    return true;
}

bool QuadraticProblem::set_q0_vector(const Vector<double> &q0) {
    if (!constraint_phase) {
        print_message("Before setting the q0 vector you should enter in the "
                      "constraint phase by either adding constraints or by "
                      "obtaining the x variable via get_main_variable()\n");
        return false;
    }
    if (q0.size() != vx->get_size()) {
        print_message("Incompatible q0 vector. q0 should be a vector of "
                      "length N with N being the size of the x variable\n");
        return false;
    }
    this->q0 = q0;
    return true;
}

double QuadraticProblem::solve_problem(Vector<double> &arg) {
    const Matrix<double> &EQ = eq.get_constraint_matrix();
    const Vector<double> &EQb = eq.get_known_values_vector();
    const Matrix<double> &GEQ = geq.get_constraint_matrix();
    const Vector<double> &GEQb = geq.get_known_values_vector();
    // pass matrices and vector in the format accepted by quadprogpp
    return solve_quadprog(Q, q0, transpose(EQ), -EQb, transpose(GEQ), -GEQb,
                          arg);
}

void QuadraticProblem::write_information() {
    std::cout << "QP: " << vx->get_size() << " variables, " <<
              eq.get_known_values_vector().size() << " equality constraints, "
              << geq.get_known_values_vector().size() <<
              " inequality constraints\n";
    pretty_print_matrix(geq.get_constraint_matrix(), "GEQ");
    print_vector("", geq.get_known_values_vector());
    pretty_print_matrix(eq.get_constraint_matrix(), "EQ");
    print_vector("", eq.get_known_values_vector());
}
