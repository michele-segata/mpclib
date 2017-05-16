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

#include "Variable.h"

void Variable::add_variable(Variable *variable) {

    // base variables cannot include subvariables
    if (base_variable) {
        std::cerr << "You can't add a variable to a base variable\n";
        return;
    }

    // get index of the previous subvariable
    int previous = variables.size() - 1;
    // this is the parent of the subvariable
    variable->parent = this;

    if (variables.size() == 0) {
        // position w.r.t. the parent variable
        variable->position = 0;
    }
    else {
        // the position is the position of the previous plus its size
        variable->position = variables[previous]->position +
                             variables[previous]->ssize;
    }
    // increment the size of this composed variable
    ssize += variable->get_size();
    // add the subvariable to the list
    variables.push_back(variable);
}

void Variable::add_variable(const Variable *variable, int n) {
    for (int i = 0; i < n; i++)
        add_variable(variable->clone());
}

int Variable::get_position() const {
    if (parent == 0)
        return position;
    else
        return parent->get_position() + position;
}

int Variable::get_size() const {
    return ssize;
}

std::string Variable::get_name() const {
    return name;
}

const Variable *Variable::get_parent() const {
    return parent;
}

Variable *Variable::get(int k) const {
    if (k >= variables.size()) {
        std::cerr << "Array index out of bound\n";
        return 0;
    }
    return variables[k];
}

std::ostream& operator<<(std::ostream &strm, const Variable &v) {
    return strm << "V.get_size() = " << v.get_size() << " V.get_position() = "
           << v.get_position();
}

std::ostream& operator<<(std::ostream &strm, const Variable *v) {
    if (v == 0)
        return strm << "NULL";
    else
        return strm << "V.get_size() = " << v->get_size() <<
               " V.get_position() = " << v->get_position();
}

void Variable::print(int level) {
    // add indentation
    for (int i = 0; i < level; i++)
        std::cout << "  ";
    // print information on this variable
    std::cout << name << " this=" << (void*)this << " parent=" <<
              (void*)parent << " position=" << get_position() << "\n";
    // print information on sub variables
    for (int i = 0; i < variables.size(); i++)
        variables[i]->print(level+1);
}

Variable *Variable::clone() const {
    // create a new variable of the same type (base_variable) with the
    // same name
    Variable *v = new Variable(base_variable, name);
    // for each subvariable the original variable has, clone it and add
    // it to the new variable
    for (int i = 0; i < variables.size(); i++)
        v->add_variable(variables[i]->clone());
    return v;
}

int Variable::get_subvariables_count() const {
    return variables.size();
}
