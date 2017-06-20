/*
 * Copyright (c) 2017 Michele Segata
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Michele Segata <segata@ccs-labs.org>
 * Description:
 *
 */
#include "MathEvaluator.h"

#include <iostream>
#include <matheval.h>
#include <string.h>

MathEvaluator::MathEvaluator(const char *expression) :
    expression(0), evaluator(0), variables(0), nVariables(0),
    values(0), assigned(0), nAssigned(0) {
	this->expression = new char[strlen(expression)];
	strcpy(this->expression, expression);
}

MathEvaluator::~MathEvaluator() {
	if (evaluator)
		evaluator_destroy(evaluator);
	if (values)
		delete [] values;
	if (assigned)
		delete [] assigned;
	if (expression)
		delete [] expression;
}

bool MathEvaluator::parseExpression() {

	evaluator = evaluator_create(expression);
	if (evaluator) {
		evaluator_get_variables(evaluator, &variables, &nVariables);
		if (nVariables > 0) {
			values = new double[nVariables];
			assigned = new bool[nVariables];
			for (int i = 0; i < nVariables; i++) {
				assigned[i] = false;
				nameToIndex[variables[i]] = i;
			}
		}
		return true;
	}
	else {
		return false;
	}

}

const char *MathEvaluator::getVariable(int i) {
	if (i >= 0 && i < nVariables)
		return variables[i];
	return 0;
}

int MathEvaluator::getVariablesCount() {
	return nVariables;
}

void MathEvaluator::printVariables(std::ostream &logFile) {
	if (evaluator) {
		int i;
		logFile << "expression: " << expression << std::endl;
		for (i = 0; i < nVariables; i++)
			logFile << "\tvariable " << i << ": " << variables[i] << std::endl;
	}
}

bool MathEvaluator::assignVariable(int i, double value) {
	if (i < 0 || i >= nVariables)
		return false;
	values[i] = value;
	if (!assigned[i]) {
		assigned[i] = true;
		nAssigned++;
	}
	return true;
}
bool MathEvaluator::assignVariable(const char *name, double value) {
	NMap::iterator v = nameToIndex.find(name);
	if (v == nameToIndex.end())
		return false;
	return assignVariable(v->second, value);
}
bool MathEvaluator::assignVariable(std::string name, double value) {
	return assignVariable(name.c_str(), value);
}
bool MathEvaluator::evaluateExpression(double &result) {

	if (!evaluator)
		return false;

	if (nAssigned != nVariables)
		return false;

	result = evaluator_evaluate(evaluator, nVariables, variables, values);
	return true;

}
bool MathEvaluator::evaluateExpression(int &result) {
	double dv;
	bool res = evaluateExpression(dv);
	result = dv;
	return res;
}
