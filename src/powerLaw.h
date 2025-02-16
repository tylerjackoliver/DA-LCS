/*
 * Created by Jack Tyler on 28/06/2021.
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

#ifndef ACROBAT_POWERLAW_H
#define ACROBAT_POWERLAW_H

#include <eigen3/Eigen/Core>
#include <stdexcept>
#include "concepts.h"
#include "linearAlgebra.h"

/* @brief Normalize a given vector.
 * @params[in] vec: Vector containing an arithmetic type.
 * @returns vector containing the normalization of the elements in vec.
 */
template <class T>
requires arithmetic<T>
std::vector<T> normalizeVector(const std::vector<T>& vec)
{
    std::vector<T> result = vec;
    T normalisingFactorSqr = 0.0, normalisingFactor;
    for (size_t i = 0; i < vec.size(); ++i)
    {
        normalisingFactorSqr += vec[i] * vec[i];
    }
    normalisingFactor = sqrt( normalisingFactorSqr );
    for (size_t i = 0; i < vec.size(); ++i)
    {
        result[i] = vec[i] / normalisingFactor;
    }
    return result;
}

/* @brief Perform power law iteration on a given matrix.
 * @param[in] mat: Matrix to determine the leading eigenvector of.
 * @param[in] vec: Approximation to the dominant eigenvector.
 * @param[in] tolerance: Optional. Iteration stops when error in consecutive steps is below this threshold.
 * @param[in] maxIterations: Optional. Iteration stops when this number of iterations is exceeded.
 */
template <class matrix, class vector>
requires has_size<matrix> and has_size<vector> and is_indexable<vector>
vector iterateToConvergence(const matrix& mat, vector& vec, double tolerance = 1e-012, int maxIterations = 10000)
{
    vector old = vec;
    double norm = std::numeric_limits<double>::max();
    int numIterations = 0;
    vector newVec;
    while (norm >= tolerance && numIterations < maxIterations) {
        norm = 0.0;
        /* Iterate & normalize. */
        newVec = normalizeVector(mat * old);                                                    // Via temporary
        vector diff(newVec.size());                                                             // Difference vector
        for (size_t i = 0; i < newVec.size(); ++i) diff[i] = newVec[i] - old[i];                // Compute diff vector
        norm = 0.0;
        /* Get relative error norm: compute the L2 norm of each order in diff and newVec, and determine the relative
         * error corresponding to each order, such that the overall norm is the maximum relative error norm.
         */
        std::vector<double> all_orders( diff[0].orderNorm().size() );
        for (auto this_element = 0; this_element < diff.size(); ++this_element){
            std::vector<double> order_norm = diff[this_element].orderNorm();
            std::vector<double> newVec_norm = newVec[this_element].orderNorm();
            for (size_t order = 0; order < order_norm.size(); ++order){
                all_orders[order] = std::max( all_orders[order], order_norm[order] / newVec_norm[order] );
            }
        }
        for (size_t i = 1; i < all_orders.size(); ++i){ // i = 1 ignores the constant part
            norm = std::max(norm, all_orders[i]);
        }
        numIterations++;
        old = newVec;
    }
    return newVec;
}

#endif //ACROBAT_POWERLAW_H
