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
#ifndef MATHEVALUATOR_H_
#define MATHEVALUATOR_H_

#include <matheval.h>
#include <iostream>
#include <map>

class MathEvaluator {
private:
	//math expression
	char   *expression;
	//evaluator instance
	void   *evaluator;
	//variables in the expression
	char  **variables;
	//number of variables in the expression
	int     nVariables;
	//values assigned to variables
	double *values;
	//tells whether values have been assigned or not
	bool   *assigned;
	//count of assigned variables
	int     nAssigned;
	//map from name to index
	typedef std::map<const char*, int> NMap;
	NMap nameToIndex;

public:
	MathEvaluator(const char *expression);
	virtual ~MathEvaluator();

	/**
	 * tries to parse the expression
	 */
	bool parseExpression();
	/**
	 * returns the name of the i-th variable, or NULL if variable i does not
	 * exists
	 */
	const char *getVariable(int i);
	/**
	 * returns the number of variables in the expression
	 */
	int getVariablesCount();
	/**
	 * logs the variables of the expression
	 */
	void printVariables(std::ostream &logFile = std::cout);
	/**
	 * sets the value of a variable
	 */
	bool assignVariable(int i, double value);
	bool assignVariable(const char *name, double value);
	bool assignVariable(std::string name, double value);
	/**
	 * evaluates the result of the mathematical expression.
	 *
	 */
	bool evaluateExpression(double &result);
	bool evaluateExpression(int &result);

};

#endif /* MATHEVALUATOR_H_ */
