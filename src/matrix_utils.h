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

#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include "Array.hh"
#include <vector>
#ifdef ENABLE_DISCRETIZATION
#include <gsl/gsl_matrix.h>
#endif

/**
 * Copies the values of a matrix over a portion of another
 * @param dst destination matrix
 * @param src source matrix
 * @param row row index of the destination matrix
 * @param col column index of the destination matrix
 */
void set_submatrix(Matrix<double> &dst, const Matrix<double> &src, int row,
                   int col);

/**
 * Copies the values of a vector over a portion of another
 * @param dst destination vector
 * @param src source vector
 * @param pos index of the destination vector
 */
void set_subvector(Vector<double> &dst, const Vector<double> &src, int pos);

/**
 * Generates a diagonal matrix by using a vector as the diagonal
 * @param src vector to use as the diagonal
 * @return the diagonal matrix
 */
Matrix<double> diag(const Vector<double> &src);

/**
 * Returns a diagonal square matrix
 * @param n square matrix size
 * @param v value to put on the diagonal
 * @return the diagonal matrix
 */
Matrix<double> diag(int n, double v = 1);

/**
 * Returns an identity matrix
 * @param n square matrix size
 * @return the identity matrix
 */
Matrix<double> identity(int n);

/**
 * Given a list of comma separated values, returns a vector
 * @param v string of comma separated values
 * @param length number of elements in the string
 * @return the vector including the parsed values
 */
Vector<double> get_vector(const std::string &v, int length);

/**
 * Returns a vector of all zeros
 * @param length length of the vector
 * @return a zero vector of the specified length
 */
Vector<double> zero_vector(int length);

/**
 * Given a string of comma separated values, returns a matrix of the parsed
 * values
 * @param m string of comma separated values. the number of elements must be
 * nrow * ncol
 * @param nrow number of rows
 * @param ncol number of columns
 * @return the matrix of parsed values
 */
Matrix<double> get_matrix(const std::string &m, int nrow, int ncol);

/**
 * Matrix transpose
 * @param m the matrix to be transposed
 * @return the transposed matrix
 */
Matrix<double> transpose(const Matrix<double> &m);

/**
 * Given a matrix of a particular size, changes the size of the matrix
 * without touching the existing values. New matrix values are initialized
 * with zeros
 * @param m matrix to be enlarged in place
 * @param r new number of rows
 * @param c new number of columns
 */
void extend_matrix(Matrix<double> &m, int r, int c);

/**
 * Given a vector of a particular size, changes the size of the vector
 * without touching the existing values. New vector values are initialized
 * with zeros
 * @param v vector to be enlarged in place
 * @param n new size
 */
void extend_vector(Vector<double> &v, int n);

/**
 * Merges two matrices together by appending the rows of the second matrix to
 * the rows of the first
 * @param a first matrix
 * @param b second matrix
 * @return a matrix where a and b are merged by rows
 */
Matrix<double> bind_matrices(const Matrix<double> &a, const Matrix<double> &b);

/**
 * Merges two vectors together
 * @param a first vector
 * @param b second vector
 * @return a vector composed by a and b
 */
Vector<double> bind_vectors(const Vector<double> &a, const Vector<double> &b);

/**
 * Performs the product of two matrices
 * @param a first matrix
 * @param b second matrix
 * @return the matrix product between a and b if the matrices are compatible,
 * an empty matrix otherwise
 */
Matrix<double> multiply(const Matrix<double> &a, const Matrix<double> &b);

/**
 * Performs the product of a matrix and a vector
 * @param a matrix
 * @param b vector
 * @return the product between a and b if the matrices are compatible,
 * an empty matrix otherwise
 */
Matrix<double> multiply(const Matrix<double> &a, const Vector<double> &b);

/**
 * Performs the product of a matrix and a vector
 * @param a matrix
 * @param b vector
 * @return the product between a and b if the matrices are compatible,
 * an empty matrix otherwise
 */
Matrix<double> multiply(const Matrix<double> &a, const double *b);

/**
 * Given a matrix, the function computes the powers of such a matrix and
 * stores it in a vector of matrices, where the index is the exponent. The
 * results is thus a vector v[i] = A^i for i = 0 ... n
 * @param a input matrix
 * @param n maximum exponent
 * @return the vector of matrices A^i
 */
std::vector<Matrix<double> > get_powers(const Matrix<double> &a, int n);

/**
 * Pretty prints a matrix to stdout
 * @param a the matrix to print
 * @param name optional name to print
 */
void pretty_print_matrix(const Matrix<double_t> &a, const std::string &name="");

/**
 * Copies a vector inside another. If the sizes don't match, then the source
 * vector is repeated multiple time in the destination vector. For example,
 * if the source vector is [1, 2, 3] and the destination vector has length 8,
 * then the destination will be filled with [1, 2, 3, 1, 2, 3, 1, 2]
 * @param dst destination vector
 * @param src source vector
 */
void copy_vector(Vector<double>& dst, const Vector<double> &src);

/**
 * Copies a vector inside another. The size of the destination vector must be
 * larger or equal than the one of the source vector
 * @param dst destination vector
 * @param src source vector
 */
void copy_vector(double *dst, const Vector<double> &src);

/**
 * Copies a vector inside another. The size of the destination vector must be
 * larger or equal than the one of the source vector
 * @param dst destination vector
 * @param src source vector
 */
void copy_vector(double *dst, const Matrix<double> &src);

/**
 * Returns a portion of a vector
 * @param src source vector
 * @param from start index
 * @param length elements to copy
 * @return the subvector of elements going from src[from] to src[from+length-1]
 */
Vector<double> subvector(const Vector<double> &src, int from, int length);

#ifdef ENABLE_DISCRETIZATION

/**
 * Converts a gsl matrix into a Matrix<double> object
 * @param m pointer to gsl_matrix
 * @return converted matrix
 */
Matrix<double> from_gsl(const gsl_matrix *m);

/**
 * Converts a Matrix<double> object into a gsl matrix. This function allocates
 * the memory for the gsl matrix and it is a duty of the user to free the
 * memory via the gsl_matrix_free() method
 * @param m matrix object
 * @return converted matrix
 */
gsl_matrix *to_gsl(const Matrix<double> &m);

#endif

#endif //MATRIX_UTILS_H
