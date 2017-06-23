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
#include "ConfigLoader.h"
#include <string>

ConfigLoader::~ConfigLoader() {}

bool ConfigLoader::getVariableName(std::string value, std::string &variable) {
    if (value.substr(0, 2) == "${" && value.substr(value.length() - 1, 1) == "}") {
        //extract variable name
        variable = value.substr(2, value.length()-3);
        return true;
    }
    return false;
}

bool ConfigLoader::parseConfig() {
    try {
        config.readFile(configFile.c_str());
        //we want automatic conversion of int to doubles
        config.setAutoConvert(true);
    }
    catch (const libconfig::ParseException &pe) {
        std::cerr << "Unable to parse " << pe.getFile() << ". Error \"" << pe.getError() << "\" at line " << pe.getLine() << std::endl;
        return false;
    }
    catch (const libconfig::FileIOException &ioe) {
        std::cerr << "Unable to parse config file. Reason: " << ioe.what() << std::endl;
        return false;
    }
    return true;
}

bool ConfigLoader::generateMatrices(const vector<double> &vA,
                                    const vector<double> &vB,
                                    const vector<double> &vC1,
                                    const vector<double> &vC2,
                                    const vector<double> &ref,
                                    const vector<double> &x0,
                                    const vector<double> &u0) {
    int n = sqrt(vA.size());
    if (n * n != vA.size()) {
        logFile << "Matrix A is not square\n";
        return false;
    }
    A.resize(n, n);
    for (int r = 0; r < n; r++)
        for (int c  = 0; c < n; c++)
            A[r][c] = vA[r * n + c];

    int p = vB.size() / n;
    if (vB.size() != n * p) {
        logFile << "Matrix B must have a number of elements which is a "
                   "multiple of n=" << n << "\n";
        return false;
    }
    B.resize(n, p);
    for (int r = 0; r < n; r++)
        for (int c  = 0; c < p; c++)
            B[r][c] = vB[r * p + c];

    int q = vC1.size() / n;
    if (vC1.size() != n * q) {
        logFile << "Matrix C1 must have a number of elements which is a "
                   "multiple of n=" << n << "\n";
        return false;
    }
    C1.resize(q, n);
    for (int r = 0; r < q; r++)
        for (int c  = 0; c < n; c++)
            C1[r][c] = vC1[r * n + c];

    int q2 = vC2.size() / n;
    if (vC2.size() != n * q2) {
        logFile << "Matrix C2 must have a number of elements which is a "
                   "multiple of n=" << n << "\n";
        return false;
    }
    C2.resize(q2, n);
    for (int r = 0; r < q2; r++)
        for (int c  = 0; c < n; c++)
            C2[r][c] = vC2[r * n + c];

    if (ref.size() == 0 || ref.size() % q != 0) {
        logFile << "The reference vector size must be a multiple of q=" << q
                << "\n";
        return false;
    }
    ref_vector.resize(ref.size());
    for (int i = 0; i < ref.size(); i++)
        ref_vector[i] = ref[i];

    if (x0.size() != n) {
        logFile << "The size of the x0 vector must be " << n << "\n";
        return false;
    }
    init_x.resize(x0.size());
    for (int i = 0; i < x0.size(); i++)
        init_x[i] = x0[i];

    if (u0.size() != p) {
        logFile << "The size of the u0 vector must be " << p << "\n";
        return false;
    }
    init_u.resize(u0.size());
    for (int i = 0; i < u0.size(); i++)
        init_u[i] = u0[i];

    return true;

}

Matrix<double> ConfigLoader::get_A() {
    return A;
}

Matrix<double> ConfigLoader::get_B() {
    return B;
}

Matrix<double> ConfigLoader::get_C1() {
    return C1;
}

Matrix<double> ConfigLoader::get_C2() {
    return C2;
}

Vector<double> ConfigLoader::get_ref_vector() {
    return ref_vector;
}

Vector<double> ConfigLoader::get_init_x() {
    return init_x;
}

Vector<double> ConfigLoader::get_init_u() {
    return init_u;
}

bool ConfigLoader::loadConfiguration() {

    vector<double> A, B, C1, C2, r, x0, u0;

    try {
        //list of simulations in config file
        Setting &cfg = config.getRoot();

        if (!recurseArray(cfg, "A", A)) {
            logFile << "missing A matrix in config file" << std::endl;
            return false;
        }
        if (!recurseArray(cfg, "B", B)) {
            logFile << "missing B matrix in config file" << std::endl;
            return false;
        }
        if (!recurseArray(cfg, "C1", C1)) {
            logFile << "missing C1 matrix in config file" << std::endl;
            return false;
        }
        if (!recurseArray(cfg, "C2", C2)) {
            logFile << "missing C2 matrix in config file" << std::endl;
            return false;
        }
        if (!recurseArray(cfg, "r", r)) {
            logFile << "missing reference vector r in config file" << std::endl;
            return false;
        }
        if (!recurseArray(cfg, "x0", x0)) {
            logFile << "missing initial state x0 in config file" << std::endl;
            return false;
        }
        if (!recurseArray(cfg, "u0", u0)) {
            logFile << "missing initial state u0 in config file" << std::endl;
            return false;
        }

        generateMatrices(A, B, C1, C2, r, x0, u0);

    }
    catch (SettingNotFoundException &snfe) {
        logFile << "No simulations found in " << configFile << std::endl;
        return false;
    }

    return true;
}

bool ConfigLoader::recurseArray(Setting &setting, const char *variable, std::vector<double> &array) {

    double value;
    try {
        //directly search for the array
        Setting &v = config.lookup(setting.getPath() + std::string(".") + std::string(variable));
        if (v.getType() == Setting::TypeArray) {
            if (v.getLength() != 0) {
                for (int i = 0; i < v.getLength(); i++) {
                    Setting &el = v[i];
                    if (el.getType() == Setting::TypeString)
                        recurseValue(el, value);
                    else
                        value = el;
                    array.push_back(value);
                }
                return true;
            }
            else {
                return false;
            }
        }
    }
    catch (SettingNotFoundException &snfe) {
        if (!setting.isRoot())
            return recurseArray(setting.getParent(), variable, array);
    }

    return false;

}
