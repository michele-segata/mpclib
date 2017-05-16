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

#include "Constraint.h"

bool Constraint::set_variable(const Var &x) {
    if (initialized || var_set)
        return false;
    this->x = x;
    var_set = true;
    // init constraint matrix
    A.resize(0, x.get_size());
    return true;
}

bool Constraint::set_constraint_variable(const Var &v,
                                         const Matrix<double> &a) {
    // sanity checks
    if (!var_set) {
        std::cerr << "Please set the main variable before adding "
                "constraint values\n";
        return false;
    }
    if (v.get_size() != a.ncols()) {
        std::cerr << "The number of columns of A matrix (" << a.ncols() <<
                  ") does not match the size of the variable (" <<
                  v.get_size() << ")\n";
        return false;
    }
    if (v.get_size() > x.get_size()) {
        std::cerr << "The size of the given variable (" << v.get_size()
                  << ") is larger than the size of the x variable ("
                  << x.get_size() << ")\n";
        return false;
    }

    if (!initialized) {
        // the size of the A matrix is set to the size of the whole x
        // variable (columns) and the size of the constraint variable v
        // (rows)
        A.resize(0, a.nrows(), x.get_size());
        initialized = true;
    } else {
        if (a.nrows() != A.nrows()) {
            std::cerr << "The number of rows of the given matrix ("
                      << a.nrows() << ") does not match the number of "
                              "rows of the A matrix (" << A.nrows() << ")\n";
            return false;
        }
    }

    set_submatrix(A, a, 0, v.get_position());
    return true;

}

bool Constraint::set_constraint_variable(const Var &v, double a) {
    Matrix<double> ma;
    ma.resize(a, 1, 1);
    return set_constraint_variable(v, ma);
}

bool Constraint::set_known_term(const Vector<double> &t) {
    // sanity check
    if (!var_set) {
        std::cerr << "Please set the main variable before adding "
                "constraint values\n";
        return false;
    }
    if (!initialized) {
        A.resize(0, t.size(), x.get_size());
    } else {
        if (t.size() != A.nrows()) {
            std::cerr << "The length of the known term vector ("
                      << t.size() << ") does not match the number of rows"
                      " of the A matrix (" << A.nrows() << ")\n";
            return false;
        }
    }

    b = t;
    return true;

}

bool Constraint::set_known_term(double t) {
    Vector<double> vt;
    vt.resize(t, 1);
    return set_known_term(vt);
}

const Matrix<double> & Constraint::get_constraint_matrix() const {
    return A;
}

const Vector<double> & Constraint::get_known_values_vector() const {
    return b;
}

bool Constraint::add_constraint(Constraint &a) {

    // sanity check
    if (a.x.get_size() != x.get_size()) {
        std::cerr << "Incompatible variable sizes when adding constraint\n";
        return false;
    }
    if (position != -1) {
        std::cerr << "WARNING: you are trying to add a constraint to a "
                  "constraint which has already been copied into the set "
                  "of constraints of the optimization problem. This can "
                  "have consequences if you decide to update the "
                  "constraint later on the optimization problem.";
    }

    if (!initialized) {
        A = a.A;
        b = a.b;
        // this is the first constraint being inserted
        a.position = 0;
        initialized = true;
    }
    else {
        // this constraint is appended at the end of the already existing
        // A matrix
        a.position = A.nrows();
        A = bind_matrices(A, a.A);
        b = bind_vectors(b, a.b);
    }
    return true;
}

bool Constraint::update_constraint(const Constraint &c, bool update_matrix,
                                   bool update_vector) {
    // has this constraint been set anywhere?
    if (c.position == -1)
        return false;
    // does the constraint fit?
    if (c.position + c.A.nrows() > A.nrows())
        return false;
    // is it compatible?
    if (c.A.ncols() != A.ncols())
        return false;
    if (update_matrix)
        set_submatrix(A, c.A, c.position, 0);
    if (update_vector)
        set_subvector(b, c.b, c.position);
    return true;
}

void Constraint::invert_constraint() {
    A = -A;
    b = -b;
}

Var Constraint::get_variable() const {
    return x;
}

void Constraint::set_position(int position) {
    this->position = position;
}
