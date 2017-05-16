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

#ifndef VARIABLE_H
#define VARIABLE_H

#include <vector>
#include <iostream>

/**
 * This is a helper class that permits to create a hierarchical set of
 * variables.
 *
 * Imagine your macro variable x is the position of a point in space and time. x
 * will be composed by a vector of variables x_k, where k indicates the time
 * sample. In turn, each x_k is composed by a vector of coordinates x_k_i, for
 * example three dimensional coordinates (x y z coordinates). Representing this
 * is not straightforward, so this class permits to define such a structure and
 * to easily obtain indexes of the single variables within the hierarchical
 * structure.
 *
 * Example:
 * x = [x_0, x_1]
 * x_0 = [x_0_0, x_0_1, x_0_2]
 * x_1 = [x_1_0, x_1_1, x_1_2]
 *
 * The full state will be [x_0_0, x_0_1, x_0_2, x_1_0, x_1_1, x_1_2]. This class
 * permits, starting from x, to obtain any index in the hierarchy. For example
 * x->get(0) refers to the block [x_0_0, x_0_1, x_0_2] within x.
 * x->get(0)->get_position() will thus return 0. x->get(1)->get_position() will
 * instead return 3. x->get(0)->get(1)->get_position() returns 1, while
 * x->get(0)->get(2)->get_position() returns 5.
 */
class Variable {

protected:
    /* indicates whether this variable is a base variable (real variable) or a
     * composed one. In the documentation example x_0_0 is a base variable,
     * while x and x_0 are not
     */
    bool base_variable;
    /* size of this variable. For a base variable, the size is 1, while for a
     * composed variable, the size is the sum of all the variables down in the
     * hierarchy
     */
    int ssize;
    /* position of the variable within the parent variable. this is set by
     * the parent variable when added with add_variable()
     */
    int position;
    /* parent variable used to recursively compute this element position. set
     * by the parent when using add_variable()
     */
    const Variable *parent;
    // name of the variable, if needed
    std::string name;
    // vector of subvariables
    std::vector<Variable *> variables;

public:

    /**
     * Constructor. Creates a new variable
     *
     * @param base_variable if true, this is a base (real) variable with size
     * one. if false, this will be a composed variable
     * @param name name of the variable
     */
    Variable(bool base_variable, const std::string &name) {
        this->base_variable = base_variable;
        this->name = name;
        // parent is set by the parent
        parent = 0;
        // the same for position
        position = 0;
        if (base_variable)
            ssize = 1;
        else
            ssize = 0;
    }

    /**
     * Destructor. Frees all the memory allocated for subvariables
     */
    ~Variable() {
        for (int i = 0; i < variables.size(); i++)
            delete variables[i];
    }

    /**
     * Adds a subvariable to this composed variable
     * @param variable the variable to add
     */
    void add_variable(Variable *variable);

    /**
     * Adds a set of subvariables to this composed variable.
     * @param variable the variable to clone multiple times to be added to
     * this composed variable. The given variable is cloned, and the class
     * doesn't take ownership of the pointer, so it is a duty of the user to
     * free the memory allocated for the variable
     * @param n number of copies to add
     */
    void add_variable(const Variable *variable, int n);

    /**
     * Returns the index of the variable within the whole hierarchy
     * @return index position
     */
    int get_position() const;

    /**
     * Returns the size of this variable. This take into account all included
     * subvariables
     * @return the size
     */
    int get_size() const;

    /**
     * Returns the name of the variable
     * @return the name
     */
    std::string get_name() const;

    /**
     * Returns the parent of the variable
     * @return the parent
     */
    const Variable *get_parent() const;

    /**
     * Returns a subvariable of this composed variable
     * @param k variable index
     * @return the subvariable of index k
     */
    Variable *get(int k) const;

    /**
     * Pretty prints this variable and all subvariables
     * @param strm the output stream on which print
     * @param v the variable
     * @return the modified output stream
     */
    friend std::ostream& operator<<(std::ostream &strm, const Variable &v);

    /**
     * Pretty prints this variable and all subvariables
     * @param strm the output stream on which print
     * @param v the variable
     * @return the modified output stream
     */
    friend std::ostream& operator<<(std::ostream &strm, const Variable *v);

    /**
     * Pretty prints the whole hierarchy structure of a variable
     * @param level subvariable level used for indenting sub variables
     */
    void print(int level = 0);

    /**
     * Clones a variable. This is useful when iteratively adding subvariables
     * which have all the same structure. Instead of creating a new one,
     * create a base one and the iteratively call add_variable() passing a
     * clone of the base variable. The function also recursively clones
     * subvariables
     * @return a pointer to the cloned instance
     */
    Variable *clone() const;

    /**
     * Returns the number of subvariables in this variable
     * @return the number of subvariables
     */
    int get_subvariables_count() const;


};

/**
 * Wrapper class for a pointer to a variable. This class is useful as it
 * overloads the [] operator, so it is easy to set constraints on
 * subvariables by referring to the variable in a math-like syntax, e.g.,
 * x[0][3].
 */
class Var{
private:
    const Variable *v;
public:
    Var(const Variable *v) {
        this->v = v;
    }
    Var operator[] (int k) const {
        return Var(v->get(k));
    }
    int get_position() const {
        return v->get_position();
    }
    int get_size() const {
        return v->get_size();
    }
    std::string get_name() const {
        return v->get_name();
    }
    Var get_parent() const {
        return Var(v->get_parent());
    }
    int get_subvariables_count() const {
        return v->get_subvariables_count();
    }
    const Variable *get_variable() const {
        return v;
    }
};

#endif //VARIABLE_H
