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

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "Array.hh"
#include "Variable.h"
#include "matrix_utils.h"

/**
 * Class for generating optimization constraints. The class keeps track of a
 * matrix A and vector b that represents either Ax = b or Ax <= b
 */
class Constraint {

private:
    // constraint matrix
    Matrix<double> A;
    // known values vector
    Vector<double> b;
    // x variable
    Var x;
    // whether A is initialized
    bool initialized;
    // whether the variable is set
    bool var_set;
    // position of this constraint within the matrix and the known values
    // vector of the main problem. this can be set when adding a constraint
    // to the main problem to track its position and edit it afterwards
    int position;

public:

    /**
     * Constructor. Takes in input the representation of the variable x for
     * setting the number of columns of the matrix A
     * @param x the x variable
     */
    Constraint(const Var &x) : x(0) {
        this->x = x;
        var_set = true;
        initialized = false;
        position = -1;
    }

    /**
     * Base constructor
     */
    Constraint() : x(0) {
        initialized = false;
        var_set = false;
        position = -1;
    }

    /**
     * Copy constructor.
     * @param c the source constraint
     */
    Constraint(const Constraint &c) : x(0) {
        x = c.get_variable();
        var_set = c.var_set;
        initialized = c.initialized;
        position = c.position;
        A = c.A;
        b = c.b;
    }

    /**
     * Sets the variable for the constraint.
     * @param x the x variable
     * @return true if the variable can be set, false if the user already
     * started adding constraint values or the variable has already been set
     */
    bool set_variable(const Var &x);

    /**
     * Adds a constraint variable to the constraint
     * @param v constraint variable
     * @param a the coefficient for the variable (a matrix of coefficients)
     * @return true if the constraint is successfully set, false otherwise
     */
    bool set_constraint_variable(const Var &v, const Matrix<double> &a);

    /**
     * Adds a constraint variable to the constraint
     * @param v constraint variable
     * @param a the coefficient for the variable
     * @return true if the constraint is successfully set, false otherwise
     */
    bool set_constraint_variable(const Var &v, double a);

    /**
     * Sets the known term of the constraint
     * @param t known term
     * @return true if the known term is successfully set, false otherwise
     */
    bool set_known_term(const Vector<double> &t);

    /**
     * Sets the known term of the constraint
     * @param t known term
     * @return true if the known term is successfully set, false otherwise
     */
    bool set_known_term(double t);

    /**
     * Returns the constraint matrix
     * @return constraint matrix A
     */
    const Matrix<double> & get_constraint_matrix() const;

    /**
     * Returns the known values vector
     * @return known values vector b
     */
    const Vector<double> & get_known_values_vector() const;

    /**
     * Adds a constraint to this constraint. The function extends the matrix
     * A and the vector b using with the matrix A and the vector b of the
     * given constraint
     * @param a the constraint to add
     * @return true if the constraint is successfully added, false if the
     * constraints are incompatible
     */
    bool add_constraint(Constraint &a);

    /**
     * Updates a sub constraint that has been included using the add_constraint
     * method. You can choose to update only the matrix, only then known values
     * vector, or both.
     * @param c the constraint
     * @param update_matrix whether to update the matrix
     * @param update_vector whether to update the known values vector
     * @return true if the update was successful, false otherwise
     */
    bool update_constraint(const Constraint &c, bool update_matrix,
                           bool update_vector);

    /**
     * Inverts the constraint. This is useful if in your problem formulation
     * you defined the inequality constraints as Ax <= b but the solver
     * accepts them in the form of Ax >= b.
     */
    void invert_constraint();

    /**
     * Returns the variable associated with the constraint
     * @return the x variable
     */
    Var get_variable() const;

    /**
     * Sets the position of the constraint in the matrix of constraints
     * @param position constraint
     */
    void set_position(int position);

};


#endif //CONSTRAINT_H
