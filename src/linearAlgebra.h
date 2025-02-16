/*
 * Created by Jack Tyler on 14/05/2021.
 * Copyright (c) 2021 University of Southampton. All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#ifndef ACROBAT_LINEARALGEBRA_H
#define ACROBAT_LINEARALGEBRA_H

#include <vector>
#include "concepts.h"
/* @brief Overloads operator* for matrix-vector multiplication with vector-of-vectors.
 * @param[in] vector<vector<T>> representing a matrix.
 * @param[in] vector<T> representing a vector.
 * @returns The result of matmul( matrix, vector ).
 * The type contained in the matrix and vector must be the same; this avoids implicit casting.
 */
template <typename T>
requires arithmetic<T>
std::vector<T> operator*(const std::vector< std::vector<T> >& matrix, const std::vector<T>& vector)
{
#ifndef RELEASE
    /* Check matrix size
     * Assume first vector is rows, each vector is columns.
     */
    if (matrix.empty() || vector.empty()) {
        throw std::runtime_error("Empty vector passed to matrix-vector multiplication.");
    }
    size_t num_cols_matrix = matrix[0].size();
    size_t num_rows_vector = vector.size();
    if (num_cols_matrix != num_rows_vector){
        throw std::runtime_error("Matrix-vector multiplication is not (m,n) x (n, 1).");
    }
    for (const auto& vecT : matrix)	{
        if ( vecT.size() != num_cols_matrix ) {
            throw std::runtime_error("Invalid matrix size in matrix-vector multiplication.");
        }
    }
#endif
    std::vector<T> result( matrix.size()  ) ;
    for (size_t i = 0; i < matrix.size(); ++i)
    {
        result[i] = 0.0;
        for (size_t j = 0; j < matrix[0].size(); ++j)
        {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    return result;
}

/* @brief Overloads the multiplication operator for vector-of-vector representations of matrices.
 * @param[in] vector<vector<T>> representing the first matrix to be multiplied.
 * @param[in] vector<vector<T>> representing the second matrix to be multiplied.
 * @returns The matrix product of the matrices.
 * Note that the type contained in the vectors must be the same in both matrices to avoid implicit casting.
 */
template <typename T>
requires arithmetic<T>
std::vector< std::vector<T> > operator*(const std::vector< std::vector<T> >& matrix1,
                                        const std::vector< std::vector<T> >& matrix2)
{
#ifndef RELEASE
    /* Check validity of inputs if compiled in non-release mode.
     * First; are either matrices empty? */
    if (matrix1.empty() || matrix2.empty()){
        throw std::runtime_error("Empty vector passed to matrix-vector multiplication.");
    }
    /* Do the dimensions agree for matrix multiplication? */
    size_t num_cols_matrix1 = matrix1[0].size();
    size_t num_rows_matrix2 = matrix2.size();
    size_t num_cols_matrix2 = matrix2[0].size();
    if (num_cols_matrix1 != num_rows_matrix2){
        throw std::runtime_error("Matrix-vector multiplication is not (m,n) x (n, o).");
    }
    /* Are the vectors representing rows consistently sized? */
    for (const auto& vecT : matrix1) {
        if (vecT.size() != num_cols_matrix1) {
            throw std::runtime_error("Invalid matrix size in matrix-matrix multiplication.");
        }
    }
    for (const auto& vecT : matrix2) {
        if ( vecT.size() != num_cols_matrix2 ) {
            throw std::runtime_error("Invalid matrix size in matrix-matrix multiplication.");
        }
    }
#endif
    std::vector< std::vector<T> > result;
    for (size_t i = 0; i < matrix1.size(); ++i)
    {
        /* Construct each row in turn & emplace_back into the vector. */
        std::vector<T> this_row( matrix2[0].size() );
        for (size_t j = 0; j < matrix2[0].size(); ++j)
        {
            /* Perform naive matmul */
            this_row[j] = 0.0;
            for (size_t k = 0; k < matrix1[0].size(); ++k)
            {
                this_row[j] += matrix1[i][k] * matrix2[k][j];
            }
        }
        /* Add row to matrix.*/
        result.emplace_back( this_row );
    }
    return result;
}

/* @brief Overloads operator* for matrix-vector multiplication with vector-of-vectors.
 * @param[in] vector<vector<T>> representing a matrix.
 * @param[in] vector<T> representing a vector.
 * @returns The result of matmul( matrix, vector ).
 * The type contained in the matrix and vector must be the same; this avoids implicit casting.
 */
template <typename T>
requires arithmetic<T>
DACE::AlgebraicVector<T> operator*(const DACE::AlgebraicVector<DACE::AlgebraicVector<T>>& matrix,
                                   const DACE::AlgebraicVector<T>& vector)
{
#ifndef RELEASE
    /* Check matrix size
     * Assume first vector is rows, each vector is columns.
     */
    if (matrix.empty() || vector.empty()) {
        throw std::runtime_error("Empty vector passed to matrix-vector multiplication.");
    }
    size_t num_cols_matrix = matrix[0].size();
    size_t num_rows_vector = vector.size();
    if (num_cols_matrix != num_rows_vector){
        throw std::runtime_error("Matrix-vector multiplication is not (m,n) x (n, 1).");
    }
    for (const auto& vecT : matrix)	{
        if ( vecT.size() != num_cols_matrix ) {
            throw std::runtime_error("Invalid matrix size in matrix-vector multiplication.");
        }
    }
#endif
    DACE::AlgebraicVector<T> result( matrix.size()  ) ;
    for (size_t i = 0; i < matrix.size(); ++i)
    {
        result[i] = 0.0;
        for (size_t j = 0; j < matrix[0].size(); ++j)
        {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    return result;
}

/* @brief Overloads the multiplication operator for vector-of-vector representations of matrices.
 * @param[in] vector<vector<T>> representing the first matrix to be multiplied.
 * @param[in] vector<vector<T>> representing the second matrix to be multiplied.
 * @returns The matrix product of the matrices.
 * Note that the type contained in the vectors must be the same in both matrices to avoid implicit casting.
 */
template <typename T>
requires arithmetic<T>
DACE::AlgebraicVector< DACE::AlgebraicVector<DACE::DA> > operator*(const DACE::AlgebraicVector< DACE::AlgebraicVector<T> >& matrix1,
                                                                   const DACE::AlgebraicVector< DACE::AlgebraicVector<T> >& matrix2)
{
#ifndef RELEASE
    /* Check validity of inputs if compiled in non-release mode.
     * First; are either matrices empty? */
    if (matrix1.empty() || matrix2.empty()){
        throw std::runtime_error("Empty vector passed to matrix-vector multiplication.");
    }
    /* Do the dimensions agree for matrix multiplication? */
    size_t num_cols_matrix1 = matrix1[0].size();
    size_t num_rows_matrix2 = matrix2.size();
    size_t num_cols_matrix2 = matrix2[0].size();
    if (num_cols_matrix1 != num_rows_matrix2){
        throw std::runtime_error("Matrix-vector multiplication is not (m,n) x (n, o).");
    }
    /* Are the vectors representing rows consistently sized? */
    for (const auto& vecT : matrix1) {
        if (vecT.size() != num_cols_matrix1) {
            throw std::runtime_error("Invalid matrix size in matrix-matrix multiplication.");
        }
    }
    for (const auto& vecT : matrix2) {
        if ( vecT.size() != num_cols_matrix2 ) {
            throw std::runtime_error("Invalid matrix size in matrix-matrix multiplication.");
        }
    }
#endif
    DACE::AlgebraicVector<DACE::AlgebraicVector<T>> result;
    for (size_t i = 0; i < matrix1.size(); ++i)
    {
        /* Construct each row in turn & emplace_back into the vector. */
        DACE::AlgebraicVector<T> this_row( matrix2[0].size() );
        for (size_t j = 0; j < matrix2[0].size(); ++j)
        {
            /* Perform naive matmul */
            this_row[j] = 0.0;
            for (size_t k = 0; k < matrix1[0].size(); ++k)
            {
                this_row[j] += matrix1[i][k] * matrix2[k][j];
            }
        }
        /* Add row to matrix.*/
        result.emplace_back( this_row );
    }
    return result;
}

/* @brief Compute the cross product of two vectors of arbitrary type and dimension.

   @param[in] vec1 The first component vector.
   @param[in] vec2 The second component vector.
   @param[out] result The cross product.
*/
template <typename T>
requires arithmetic<T>
void cross3( const std::vector<T>& vec1, const std::vector<T>& vec2, std::vector<T>& result )
{
    /* Make sure vec1 and vec2 are of the same (correct) size.*/
#ifndef RELEASE
    if ( vec1.size() != 3 ) throw std::runtime_error("One of vec1 or vec2 are not of the correct size.");
    if ( vec1.size() != vec2.size() ) throw std::runtime_error("vec1 and vec2 are not the same size in cross.");
#endif
    /* Resize result to ensure correct size & memory allocated. */
    result.resize( vec1.size() );
    /* Perform cross product */
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

/* @brief Compute the cross product of two vectors of arbitrary type and dimension.

   @param[in] vec1 The first component vector.
   @param[in] vec2 The second component vector.
   @param[out] result The cross product.
*/
template <typename T>
requires arithmetic<T>
DACE::AlgebraicVector<T> cross3( const DACE::AlgebraicVector<T>& vec1, const DACE::AlgebraicVector<T>& vec2 )
{
    /* Make sure vec1 and vec2 are of the same (correct) size.*/
#ifndef RELEASE
    if ( vec1.size() != 3 ) throw std::runtime_error("One of vec1 or vec2 are not of the correct size.");
    if ( vec1.size() != vec2.size() ) throw std::runtime_error("vec1 and vec2 are not the same size in cross.");
#endif
    /* Resize result to ensure correct size & memory allocated. */
    DACE::AlgebraicVector<T> result(3);
    /* Perform cross product */
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    return result;
}

/* @brief Compute the dot product of two arbitrarily-sized vectors.
   @param[in] vec1 First component vector.
   @param[in] vec2 Second component vector.
*/
template <typename T>
requires arithmetic<T>
T dotProduct( const std::vector<T>& vec1, const std::vector<T>& vec2 )
{
#ifndef RELEASE
    if ( vec1.empty() || vec2.empty() ){
        throw std::runtime_error("You passed an empty vector into dotProduct; result undefined!");
    }
    if ( vec1.size() != vec2.size() ){
        throw std::runtime_error("One of vec1 or vec2 are not the same size in dotProduct.");
    }
#endif
    T result = 0.0;
    for (size_t idx = 0; idx < vec1.size(); idx++)
    {
        result += vec1[idx] * vec2[idx];
    }
    return result;
}

/* @brief Compute the Euclidean vector norm (length) for a given vector.
   @param[in] Desired vector to compute the norm for.
   @return Euclidean norm.
*/
template <typename T>
requires arithmetic<T>
T vectorNorm( const std::vector<T>& vector )
{
#ifndef RELEASE
    if (vector.empty()) throw std::runtime_error("You passed an empty vector to vectorNorm.");
#endif
    T result = 0.0;
    for (size_t idx = 0; idx < vector.size(); ++idx)
    {
        result += vector[idx] * vector[idx];
    }
    return sqrt( result );
}

/* @brief Compute the normalized version of an input vector.
   @param[in] Vector to compute the norm in.
   @param[out] Normalised vector
*/
template <typename T>
requires arithmetic<T>
void normalizeVector( std::vector<T> vec1, std::vector<T>& unitVector )
{
#ifndef RELEASE
    /* Check inputs. */
    if ( vec1.empty() ) {
        throw std::runtime_error("You passed an empty vector into normalizeVector!");
    }
#endif
    /* Ensure unit vector is the correct size. */
    unitVector.resize( vec1.size() );
    /* Get vector norm */
    T norm = vectorNorm( vec1 );
    for (size_t idx = 0; idx < vec1.size(); ++idx)
    {
        unitVector[idx] = vec1[idx] / norm;
    }
}

/* @brief Compute the normalized version of an input vector.
   @param[inout] Vector to compute the norm of.
*/
template <typename T>
requires arithmetic<T>
void normalizeVector( std::vector<T>& vec1 )
{
#ifndef RELEASE
    /* Check inputs. */
    if ( vec1.empty() ) {
        throw std::runtime_error("You passed an empty vector into normalizeVector!");
    }
#endif
    std::vector<T> unitVector( vec1.size() );
    /* Ensure unit vector is the correct size. */
    unitVector.resize( vec1.size() );
    /* Get vector norm */
    T norm = vectorNorm( vec1 );
    for (size_t idx = 0; idx < vec1.size(); ++idx)
    {
        unitVector[idx] = vec1[idx] / norm;
    }
    vec1 = unitVector;
}

/* @brief Check if a given matrix 'in' is square.
 * @param[in] in: Matrix to check.
 * @returns bool: True/False if square.
 */
template <typename T> requires is_indexable<T>
bool is_square(const T& in){
    int num_cols = in[0].size();
    int num_rows = in.size();
    if (num_cols != num_rows) throw std::runtime_error("Matrix is not square!");
    for (auto i = 0; i < in.size(); ++i) {
        if (in[i].size() != num_cols) throw std::runtime_error("Matrix is unevenly sized!");
    }
}

/* @brief Transpose a matrix comprised of AlgebraicVectors-of-AlgebraicVectors.
 * @param[in] AlgebraicVector<AlgebraicVector<DA>> in: Matrix to transpose.
 * @returns Transposed matrix corresponding to the input.
 */
DACE::AlgebraicVector<DACE::AlgebraicVector<DACE::DA>> transpose(const DACE::AlgebraicVector<DACE::AlgebraicVector<DACE::DA>>& in){
#ifndef RELEASE
    if (in.empty()) throw std::runtime_error("Empty vector passed to transpose.");
    for (size_t i = 0; i < in.size(); ++i) if (in[i].empty()) throw std::runtime_error("Bad matrix passed to transpose.");
#endif
    size_t num_cols = in[0].size();
    size_t num_rows = in.size();
    DACE::AlgebraicVector< DACE::AlgebraicVector<DACE::DA> > out;
    for (size_t i = 0; i < num_cols; ++i){
        DACE::AlgebraicVector<DACE::DA> to_append(num_rows);
        out.emplace_back(to_append);
    }
    for (size_t i = 0; i < in.size(); ++i){
        for (size_t j = 0; j < in[0].size(); ++j){
            out[j][i] = in[i][j];
        }
    }
    return out;
}

/* @brief Multiply an AlgebraicVector element-wise containing arbitrary type with another type.
 * @param[in] vec: AlgebraicVector<T> to multiply element-wise.
 * @param[in] coefficient: Type<T> to multiply.
 * @returns The result of vec * coefficient.
 */
template <typename T>
requires arithmetic<T>
DACE::AlgebraicVector<T> operator*(const DACE::AlgebraicVector<T>& vec, const T& coefficient){
    DACE::AlgebraicVector<T> new_vec( vec.size() );
    for (size_t i = 0; i < new_vec.size(); ++i){
        new_vec[i] = vec[i] * coefficient;
    }
    return new_vec;
}

/* @brief Overload operator* for coefficient * AlgebraicVector.
 * @param[in] coefficient: Type<T> to multiply the AlgebraicVector by.
 * @param[in] vec: AlgebraicVector<T> to multiply.
 * @returns The result of coefficient * vec.
 */
template <typename T, typename U>
requires arithmetic<T> and arithmetic<U>
DACE::AlgebraicVector<T> operator*(const T& coefficient, const DACE::AlgebraicVector<T>& vec){
    DACE::AlgebraicVector<T> new_vec(vec.size());
    for (size_t i = 0; i < new_vec.size(); ++i){
        new_vec[i] = vec[i] * coefficient;
    }
    return new_vec;
}

/* @brief Compute the element-wise multiplication of a matrix by a given coefficient.
 * @param[in] coefficient: Type<T> to multiply by.
 * @param[in] matrix: AlgebraicVector<AlgebraicVector<T>> to multiply.
 * @returns The result of coefficient * matrix (matrix-wise.)
 */
template <typename T, typename U>
requires arithmetic<T> and arithmetic<U>
DACE::AlgebraicVector<DACE::AlgebraicVector<DACE::DA>> operator*(const T* coefficient, const DACE::AlgebraicVector<DACE::AlgebraicVector<T>>& matrix){
    DACE::AlgebraicVector<DACE::AlgebraicVector<DACE::DA>> out = matrix;
    for (size_t i = 0; i < out.size(); ++i){
        for (size_t j = 0; j < out[0].size(); ++j){
            out[i][j] = matrix[i][j] * coefficient;
        }
    }
    return out;
}

/* @brief Perform the element-wise multiplication of a matrix by a coefficient.
 * @param[in] matrix: AlgebraicVector<AlgebraicVector<T>> to multiply.
 * @param[in] coefficient: Type<T> to multiply matrix by.
 * @returns The result of matrix * coefficient.
 */
template <typename T, typename U>
requires arithmetic<T> and arithmetic<U>
DACE::AlgebraicVector<DACE::AlgebraicVector<DACE::DA>> operator*(const DACE::AlgebraicVector<DACE::AlgebraicVector<T>>& matrix, const T* coefficient){
    DACE::AlgebraicVector<DACE::AlgebraicVector<DACE::DA>> out = matrix;
    for (size_t i = 0; i < out.size(); ++i){
        for (size_t j = 0; j < out[0].size(); ++j){
            out[i][j] = matrix[i][j] * coefficient;
        }
    }
    return out;
}

/* @brief Overload operator/ to perform element-wise division of a vector by a coefficient.
 * @param[in] vec: Element to perform the element-wise division on.
 * @param[in] coefficient: Object to divide by. Same type as those contained in vec.
 * @returns AlgebraicVector<T> corresponding to result of vec / coefficient.
 */
template <typename T>
requires arithmetic<T>
DACE::AlgebraicVector<T> operator/(const DACE::AlgebraicVector<T>& vec, const T& coefficient){
    DACE::AlgebraicVector<T> new_vec = vec;
    for (size_t i = 0; i < new_vec.size(); ++i){
        new_vec[i] /= coefficient;
    }
    return new_vec;
}

/* @brief Overload operator/ to perform element-wise division of a vector by a coefficient.
 * @param[in] coefficient: Object to divide by. Same type as those contained in vec.
 * @param[in] vec: Element to perform the element-wise division on.
 * @returns AlgebraicVector<T> corresponding to result of vec / coefficient.
 */
template <typename T>
requires arithmetic<T>
DACE::AlgebraicVector<T> operator/(const T& coefficient, const DACE::AlgebraicVector<T>& vec){
    DACE::AlgebraicVector<T> new_vec = vec;
    for (size_t i = 0; i < new_vec.size(); ++i){
        new_vec[i] /= coefficient;
    }
    return new_vec;
}

/* @brief Overload operator<< to allow printing of arbitrary vectors.
 * @param[in] ostream object.
 * @param[in] vec: Vector to print the elements of.
 * @returns Reference to the ostream object.
 */
template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T>& vec)
{
    for (size_t i = 0; i < vec.size()-1; ++i)
    {
        out << vec[i] << ", ";
    }
    return out << vec[ vec.size() - 1 ];
}

#endif //ACROBAT_LINEARALGEBRA_H
