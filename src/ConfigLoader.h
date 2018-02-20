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

#ifndef CONFIGLOADER_H_
#define CONFIGLOADER_H_

#include <libconfig.h++>
using namespace libconfig;
#include <iostream>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <string.h>
#include "MathEvaluator.h"
#include "Array.hh"
using namespace std;

/**
 * Loads simulations from configuration file
 */
class ConfigLoader {

private:
    //input file
    std::string configFile;
    //config parser
    Config config;
    //output log
    std::ostream &logFile;
    //loaded configuration
    Matrix<double> A, B, C1, C2;
    Vector<double> ref_vector, init_x, init_u;

    /**
     * Return, if any, the name of a variable inside ${}
     */
    bool getVariableName(std::string value, std::string &variable);

    /**
     * recursively search value for a variable, evaluating variables as well,
     * e.g., ${parameter}
     *
     * \param setting the setting where to start searching from. if not found
     * in here, function will search in the parents
     * \param variable name of the variable to be search
     * \param result where to store the final result
     * \return true if the value is found, false otherwise
     */
    template <typename T>
    bool recurseVariable(Setting &setting, const char *variable, T &result) {

        std::string strValue;
        std::string var;

        //search for the variable called "variable"
        if (setting.lookupValue(variable, strValue)) {
            //we found the value, is this a variable?
            if (getVariableName(strValue, var)) {
                //if it is a variable, check whether the name is the same as the
                //requested one. if so, we need to go search in the parent to
                //avoid infinite recursion
                if (var == variable)
                    return recurseVariable(setting.getParent(), var.c_str(), result);
                //otherwise search again in the current setting
                else
                    return recurseVariable(setting, var.c_str(), result);
            }
            else {
                //this is not a variable, but it is a string, so we return the value
                //here we call again lookup value because we actually don't know the
                //type of "result". here we know that strValue is a string but maybe
                //the caller was searching for an int
                return setting.lookupValue(variable, result);
            }
        }
        //we did not find the variable
        else {
            //search its value directly. note: if we're here, this is not a string
            //so we don't search for a variable
            if (setting.lookupValue(variable, result)) {
                return true;
            }
            else {
                //we need to go search for it in the parent
                if (!setting.isRoot())
                    return recurseVariable(setting.getParent(), variable, result);
                else
                    return false;
            }
        }

    }

    template <typename T>
    bool recurseValue(Setting &setting, T &result) {

        std::string strValue;
        std::string var;

        //search for the variable called "variable"
        if (setting.getType() == setting.TypeString) {
            strValue = std::string(setting.c_str());
            //we found the value as a string. invoke a math evaluator
            MathEvaluator me(strValue.c_str());
            //check if the expression is fine
            if (!me.parseExpression())
                return false;
            //is there any variable?
            if (me.getVariablesCount() == 0)
                //if not we can simply evaluate the expression
                return me.evaluateExpression(result);

            //if there are variables, then we need to find them
            for (int i = 0; i < me.getVariablesCount(); i++) {
                double variableValue;
                bool searchResult;
                if (!recurseExpression(setting, me.getVariable(i), variableValue))
                    //variable was not found in config file
                    return false;

                //variable was found, assign its value
                me.assignVariable(i, variableValue);

            }

            //we have found and assigned all the variables. evaluate the expression
            return me.evaluateExpression(result);

        }

        return false;

    }

    template <typename T>
    bool recurseExpression(Setting &setting, const char *variable, T &result) {

        std::string strValue;
        std::string var;

        //search for the variable called "variable"
        if (setting.lookupValue(variable, strValue)) {
            //we found the value as a string. invoke a math evaluator
            MathEvaluator me(strValue.c_str());
            //check if the expression is fine
            if (!me.parseExpression())
                return false;
            //is there any variable?
            if (me.getVariablesCount() == 0)
                //if not we can simply evaluate the expression
                return me.evaluateExpression(result);

            //if there are variables, then we need to find them
            for (int i = 0; i < me.getVariablesCount(); i++) {
                double variableValue;
                bool searchResult;
                //if we search for a variable with the same name as the requested
                //one, then go one level up, otherwise we'll incur in infinite
                //recursion
                if (strcmp(variable, me.getVariable(i)) == 0) {
                    if (!recurseExpression(setting.getParent(), me.getVariable(i), variableValue))
                        //variable was not found in config file
                        return false;
                }
                else {
                    //if the variable name is different, start searching at the same level
                    if (!recurseExpression(setting, me.getVariable(i), variableValue))
                        //variable was not found in config file
                        return false;
                }

                //variable was found, assign its value
                me.assignVariable(i, variableValue);

            }

            //we have found and assigned all the variables. evaluate the expression
            return me.evaluateExpression(result);

        }
        //we did not find the variable
        else {
            //search its value directly. note: if we're here, this is not a string
            //so we don't search for a variable
            if (setting.lookupValue(variable, result)) {
                return true;
            }
            else {
                //we need to go search for it in the parent
                if (!setting.isRoot())
                    return recurseExpression(setting.getParent(), variable, result);
                else
                    return false;
            }
        }

    }

    /**
     * recursively searches an array
     */
    bool recurseArray(Setting &setting, const char *variable, std::vector<double> &array);

public:

    ConfigLoader(std::string configFile, std::ostream &logFile = std::cerr) :
        configFile(configFile), logFile(logFile) {}
    ConfigLoader(const char *configFile, std::ostream &logFile = std::cerr) :
        configFile(configFile), logFile(logFile) {}
    virtual ~ConfigLoader();

    /**
     * Parses the config file
     */
    bool parseConfig();

    /**
     * Loads the configuration for all simulations
     */
    bool loadConfiguration();

    /**
     * Generates system matrices starting from loaded vectors
     * @param vA vector representing the continuous time state matrix
     * @param vB vector representing the continuous time control matrix
     * @param vC1 vector representing output matrix C1
     * @param vC2 vector representing output matrix C2
     * @param ref reference vector
     */
    bool generateMatrices(const vector<double> &vA, const vector<double> &vB,
                          const vector<double> &vC1, const vector<double> &vC2,
                          const vector<double> &ref, const vector<double> &x0,
                          const vector<double> &u0);

    /**
     * Gets the discretized state matrix
     * @return state matrix
     */
    Matrix<double> get_A();

    /**
     * Gets the discretized control matrix
     * @return control matrix
     */
    Matrix<double> get_B();

    /**
     * Gets the discretized C1 output matrix
     * @return C1 output matrix
     */
    Matrix<double> get_C1();

    /**
     * Gets the discretized C2 output matrix
     * @return C2 output matrix
     */
    Matrix<double> get_C2();

    /**
     * Gets the reference vector
     * @return reference vector
     */
    Vector<double> get_ref_vector();

    /**
     * Gets the initial state
     * @return initial state
     */
    Vector<double> get_init_x();

    /**
     * Gets the initial control
     * @return initial control
     */
    Vector<double> get_init_u();

};

#endif /* CONFIGLOADER_H_ */
