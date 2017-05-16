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

#include "matrix_utils.h"

#include <sstream>

void set_submatrix(Matrix<double> &dst, const Matrix<double> &src, int row,
                   int col) {
    int i, j;
    for (i = row; i < row + src.nrows(); i++) {
        for (j = col; j < col + src.ncols(); j++) {
            dst[i][j] = src[i - row][j - col];
        }
    }
}

void set_subvector(Vector<double> &dst, const Vector<double> &src, int pos) {
    int i;
    for (i = pos; i < pos + src.size(); i++)
        dst[i] = src[i - pos];
}

Matrix<double> diag(const Vector<double> &src) {
    Matrix<double> d;
    d.resize(0, src.size(), src.size());
    int i;
    for (i = 0; i < src.size(); i++)
        d[i][i] = src[i];
    return d;
}

Matrix<double> diag(int n, double v) {
    Vector<double> iv;
    iv.resize(v, n);
    return diag(iv);
}

Matrix<double> identity(int n) {
    return diag(n, 1);
}

Vector<double> get_vector(const std::string &v, int length) {
    char ch;
    Vector<double> rv;
    rv.resize(length);
    std::istringstream iss(v);
    for (int i = 0; i < length; i++)
        iss >> rv[i] >> ch;
    return rv;
}

Vector<double> get_vector(double v, int length) {
    Vector<double> rv;
    rv.resize(v, length);
    return rv;
}

Vector<double> zero_vector(int length) {
    return get_vector(0, length);
}

Matrix<double> get_matrix(const std::string &m, int nrow, int ncol) {
    char ch;
    Matrix<double> rm;
    rm.resize(nrow, ncol);
    std::istringstream iss(m);
    for (int i = 0; i < nrow; i++)
        for (int j = 0; j < ncol; j++)
            iss >> rm[i][j] >> ch;
    return rm;
}

Matrix<double> transpose(const Matrix<double> &m) {
    Matrix<double> t;
    t.resize(m.ncols(), m.nrows());
    for (int r = 0; r < m.nrows(); r++)
        for (int c = 0; c < m.ncols(); c++)
            t[c][r] = m[r][c];
    return t;
}

void extend_matrix(Matrix<double> &m, int r, int c) {

    if (r < m.nrows() || c < m.ncols())
        return;

    Matrix<double> tmp = m;
    m.resize(0, r, c);
    set_submatrix(m, tmp, 0, 0);

}

void extend_vector(Vector<double> &v, int n) {
    if (n < v.size())
        return;

    Vector<double> tmp = v;
    v.resize(0, n);
    set_subvector(v, tmp, 0);

}

Matrix<double> bind_matrices(const Matrix<double> &a,
                             const Matrix<double> &b) {
    Matrix<double> m;
    if (a.ncols() != b.ncols()) {
        std::cerr << "Cannot bind matrices. Mismatching column number\n";
    }
    else {
        m = a;
        extend_matrix(m, a.nrows() + b.nrows(), a.ncols());
        set_submatrix(m, b, a.nrows(), 0);
    }
    return m;
}

Vector<double> bind_vectors(const Vector<double> &a, const Vector<double> &b) {
    Vector<double> v;
    v = a;
    extend_vector(v, a.size() + b.size());
    set_subvector(v, b, a.size());
    return v;
}

Matrix<double> multiply(const Matrix<double> &a, const Matrix<double> &b) {
    Matrix<double> res;
    if (a.ncols() != b.nrows()) {
        std::cerr << "Incompatible matrix sizes. Cannot multiply\n";
        return res;
    }
    res.resize(0, a.nrows(), b.ncols());
    // non optimized multiplication (O(n^3))
    for (int r = 0; r < a.nrows(); r++)
        for (int c = 0; c < b.ncols(); c++)
            for (int i = 0; i < a.ncols(); i++)
                res[r][c] += a[r][i] * b[i][c];
    return res;
}

Matrix<double> multiply(const Matrix<double> &a, const Vector<double> &b) {
    Matrix<double> v(1, b.size());
    copy_vector(v[0], b);
    return multiply(a, transpose(v));
}

Matrix<double> multiply(const Matrix<double> &a, const double *b) {
    Vector<double> v(b, a.ncols());
    return multiply(a, v);
}

std::vector<Matrix<double> > get_powers(const Matrix<double> &a, int n) {
    std::vector<Matrix<double> > res;

    if (a.nrows() != a.ncols()) {
        std::cerr << "Matrix is not squared. Can't compute powers" << "\n";
        return res;
    }
    // A^0 = I
    res.push_back(identity(a.nrows()));
    // A^1 = A
    res.push_back(a);
    for (int i = 2; i <= n; i++)
        res.push_back(multiply(res[i-1], a));
    return res;
}

void pretty_print_matrix(const Matrix<double_t> &a, const std::string &name) {
    if (name.compare("") != 0)
        std::cout << name << ": ";
    std::cout << a.nrows() << " rows x " << a.ncols() << " columns\n";
    std::vector<int> max_lens;
    for (int c = 0; c < a.ncols(); c++) {
        max_lens.push_back(0);
        for (int r = 0; r < a.nrows(); r++) {
            std::ostringstream strs;
            strs << a[r][c];
            std::string str = strs.str();
            max_lens[c] = str.size() > max_lens[c] ? str.size() :
                          max_lens[c];
        }
    }
    // save default formatting
    for (int r = 0; r < a.nrows(); r++) {
        for (int c = 0; c < a.ncols(); c++) {
            std::ostringstream strs;
            strs << a[r][c];
            std::string str = strs.str();
            for (int ws = 0; ws < max_lens[c]-str.size(); ws++)
                std::cout << " ";
            std::cout << str << ", ";
        }
        std::cout << "\n";
    }

}

void copy_vector(Vector<double>& dst, const Vector<double> &src) {
    for (int i = 0; i < dst.size(); i++)
        dst[i] = src[i % src.size()];
}

void copy_vector(double *dst, const Vector<double> &src) {
    for (int i = 0; i < src.size(); i++)
        dst[i] = src[i];
}

void copy_vector(double *dst, const Matrix<double> &src) {
    // assume the vector to be a matrix with a single column
    for (int i = 0; i < src.nrows(); i++)
        dst[i] = src[i][0];
}

Vector<double> subvector(const Vector<double> &src, int from, int length) {
    Vector<double> sub(length);
    for (int i = 0; i < length; i++)
        sub[i] = src[i+from];
    return sub;
}
