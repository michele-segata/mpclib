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

#ifndef QUADRATICPROBLEM_H
#define QUADRATICPROBLEM_H

#include "Variable.h"
#include "Constraint.h"

#include <map>

/**
 * Defines a standard interface for the definition of quadratic optimization
 * problems that will be solved using the QuadProg++ library. The form of the
 * problem is
 *
 * min 1/2 x^T Q x + q0 x
 *
 * subject to
 *
 * Ax >= b
 * Cx <= d
 * Ex = f
 *
 * where x is the variable (a column vector of variables), Q is a symmetric
 * matrix used for weighting variables during the minimization, q0 is the vector
 * of linear coefficients, A and b are the matrix and the vector of known terms
 * for inequality (>=) constraints, C and d are the matrix and the vector of
 * known terms for inequality (<=) constraints, and E and f are the matrix and
 * the vector of known terms for equality constraints.
 *
 * Supposing that the number of variables is N, the number of inequality
 * constraints (>=) is NI, the number of inequality constraints (<=) is NI',
 * and the number of equality constraints is NE, then Q is an N x N diagonal
 * matrix, A is an N x NI matrix and b a vector of length NI, B is an N x NI'
 * matrix and d a vector of length NI', and E is an N x NE matrix and f a
 * vector of length NE.
 *
 * This class is first used to create the variable x. The class provides
 * utilities to generate variables and subvariables and add them to the main
 * x variable of the optimization problem. Once this step is completed, the
 * class can be used to insert constraints. The class automatically checks
 * the compatibility of the given constraints, i.e., it checks whether given
 * matrices have the correct size w.r.t. to the variable x, so once the user
 * starts inserting constraints, it is not possible to change the problem
 * variable x anymore.
 *
 */
class QuadraticProblem {

protected:
    // variable of the problem
    Variable *vx;
    // equality constraints
    Constraint eq;
    // inequality constraints. we only need one type of constraints. <=
    // constraints are automatically converted into >= constraints
    Constraint geq;
    // Q matrix
    Matrix<double> Q;
    // q0 vector
    Vector<double> q0;
    // indicates whether the user started inserting constraints
    bool constraint_phase;
    // indicates whether debug messages should be printed
    bool debug;

    // used to keep track of the variables that the user added to the main x
    // variable
    typedef std::map<Variable *, int> VarMap;
    VarMap variables;

    /**
     * If the debug flag is enabled, prints a message to stderr
     * @param msg the message to be printed
     */
    void print_message(const std::string &msg);

    /**
     * Starts the constraint adding phase
     */
    void init_constraint_phase();

public:

    QuadraticProblem(bool debug_output = false) {
        vx = new Variable(false, "x");
        constraint_phase = false;
        debug = debug_output;
    }

    ~QuadraticProblem() {
        delete vx;
    }

    /**
     * Generates a vector variable. For example, if you have a variable x
     * which is composed by subvariables x_k, k = 0..N-1, you can call this
     * method to generate x
     * @param n number of variables in the vector
     * @param name name to be given to the variable (e.g., x)
     * @param base_name name to be given to subvariables (e.g., x_k)
     * @return the generated vector variable
     */
    Variable *vector_variable(int n, const std::string &name = "",
                              const std::string &base_name = "");

    /**
     * Generates a vector variable of variables.
     * @param n number of variables in the vector
     * @param base_variable the base variable to be used as subvariable
     * @param name name to be given to the variable
     * @return the generated vector variable
     */
    Variable *vector_variable(int n, const Variable *base_variable,
                              const std::string &name = "");

    /**
     * Adds a variable to the main variable x
     * @param v the variable to be added.
     * @return true if successfull, false if the variable cannot be added
     * because the user already started adding constraints or because it has
     * already been added
     */
    bool add_variable(Variable *v);

    /**
     * Converts a pointer to Variable to a Var instance, which is useful for
     * adding constraints by using the operator [].
     * @param v the variable
     * @return the Var instance associated to the given v variable
     */
    Var get_variable(const Variable *v);

    /**
     * Returns the x variable of the problem, which can be used, if needed,
     * to generate constraints based on the variable. Calling this function
     * will move the class to the constraint phase, disabling any possibility
     * to add additional variable to the problem
     * @return Var instance for x
     */
    Var get_main_variable();

    /**
     * Adds an equality constraint to the problem. The constraint is passed
     * by reference because the function will update the position field of
     * the constraint to keep track of its position within the whole set of
     * constraints, enabling later modification.
     * @param eq the constraint
     * @return true if successfull, false if not. A reason for failure is the
     * fact that the given constraint is not built for the x variable
     * associated with this problem
     */
    bool add_equality_constraint(Constraint &c);

    /**
     * Updates an equality constraint. You can choose to update only the
     * matrix, only then known values vector, or both.
     * @param c the constraint
     * @param update_matrix whether to update the matrix
     * @param update_vector whether to update the known values vector
     * @return true if the update was successful, false otherwise
     */
    bool update_equality_constraint(const Constraint &c, bool update_matrix,
                                    bool update_vector);

    /**
     * Adds a >= constraint. The constraint is passed by reference because the
     * function will update the position field of the constraint to keep track
     * of its position within the whole set of constraints, enabling later
     * modification.
     * @param eq the constraint
     * @return true if successfull, false if not. A reason for failure is the
     * fact that the given constraint is not built for the x variable
     * associated with this problem
     */
    bool add_geq_constraint(Constraint &c);

    /**
     * Updates a >= constraint. You can choose to update only the matrix, only
     * then known values vector, or both.
     * @param c the constraint
     * @param update_matrix whether to update the matrix
     * @param update_vector whether to update the known values vector
     * @return true if the update was successful, false otherwise
     */
    bool update_geq_constraint(const Constraint &c, bool update_matrix,
                               bool update_vector);

    /**
     * Adds a <= constraint. The constraint is passed by reference because the
     * function will update the position field of the constraint to keep track
     * of its position within the whole set of constraints, enabling later
     * modification.
     * @param eq the constraint
     * @return true if successfull, false if not. A reason for failure is the
     * fact that the given constraint is not built for the x variable
     * associated with this problem
     */
    bool add_leq_constraint(Constraint &c);

    /**
     * Updates a <= constraint. You can choose to update only the matrix, only
     * then known values vector, or both.
     * @param c the constraint
     * @param update_matrix whether to update the matrix
     * @param update_vector whether to update the known values vector
     * @return true if the update was successful, false otherwise
     */
    bool update_leq_constraint(const Constraint &c, bool update_matrix,
                               bool update_vector);
    /**
     * Sets the Q matrix of the problem
     * @param Q the matrix
     * @return true if the matrix was set, false if the operation failed. The
     * operation might fail because the class is not already in constraint
     * mode or because the size of the matrix is incompatible
     */
    bool set_Q_matrix(const Matrix<double> &Q);

    /**
     * Sets the q0 vector of the problem
     * @param q0 the vector
     * @return true if the vector was set, false if the operation failed. The
     * operation might fail because the class is not already in constraint
     * mode or because the size of the vector is incompatible
     */
    bool set_q0_vector(const Vector<double> &q0);

    /**
     * Solves the optimization problem.
     * @param arg vector of doubles where the argument that minimizes the
     * function is returned
     * @return the minimum value
     */
    double solve_problem(Vector<double> &arg);

    /**
     * Dumps some problem information on stdout
     */
    void write_information();

};


#endif //QUADRATICPROBLEM_H
