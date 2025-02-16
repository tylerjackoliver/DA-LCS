/*
 * Created by Jack Tyler on 07/06/2021.
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

#ifndef ACROBAT_EIGENVECTORS_H
#define ACROBAT_EIGENVECTORS_H

#include <dace/dace_s.h>
#include <dace/AlgebraicVector.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <cmath>
#include "linearAlgebra.h"
#include "rhs.h"
#include "integration.h"
#include <vector>
#include "powerLaw.h"
#include "er3bpMapping.h"
#include <filesystem>

/* @brief Adjust a given matrix of AlgebraicVector-of-AlgebraicVectors to a certain size.
 * @param[in] x0: The matrix to adjust.
 * @param[in] num_rows: Number of rows of the matrix.
 * @param[in] num_cols: Number of columns of the matrix.
 * @returns Matrix with given rows and columns.
 */
auto makeMatrixToSize(DACE::AlgebraicVector<DACE::AlgebraicVector<DACE::DA>>& x0,
                      int num_rows, int num_cols) -> void {
    /* Empty x0 first */
    x0.clear();
    for (auto i = 0; i < num_rows; i++){
        DACE::AlgebraicVector<DACE::DA> tmp(num_cols);
        x0.emplace_back(tmp);
    }
}

/* @brief Create a matrix of AlgebraicVector-of-AlgebraicVectors to a certain size.
 * @param[in] num_rows: Number of rows of the matrix.
 * @param[in] num_cols: Number of columns of the matrix.
 * @returns Matrix with given rows and columns.
 */
auto makeMatrix(int num_rows, int num_cols) -> DACE::AlgebraicVector<DACE::AlgebraicVector<DACE::DA>>{
    DACE::AlgebraicVector<DACE::AlgebraicVector<DACE::DA>> to_return;
    makeMatrixToSize(to_return, num_rows, num_cols);
    return to_return;
}

/* @brief Compute the flow corresponding to an initial condition over a given time.
 * @param[in] rthphi: The initial condition.
 * @param[in] t0: The initial time.
 * @param[in] tf: The final time.
 * @returns The flow map of the initial condition over the given time.
 */
auto getFlowMap(const DACE::AlgebraicVector<DACE::DA>& rthphi, double t0, double tf,
				bool save = false, int id = -1) -> DACE::AlgebraicVector<DACE::DA>
{
	/* Example flow map using orbital elements */
	using DACE::AlgebraicVector;
	using DACE::DA;
	AlgebraicVector<DA> initial_spherical = rthphi + AlgebraicVector<DA>::identity();
	AlgebraicVector<DA> initial_oe = {initial_spherical[0], PARAMS::ECC, initial_spherical[1], PARAMS::LONGTD, initial_spherical[2], 0.0};
	AlgebraicVector<DA> initial_position(6); oes_to_state_simplified(initial_oe, PARAMS::targetGM, initial_position);
    AlgebraicVector<DA> er3bp_state = secondary_to_er3bp(initial_position, t0);
	integration<DA, double> integrator(er3bp_state, t0, tf);
	integrator.setTol(1e-013);
	if (id >= 0) {
		integrator.setId(id);
	}
	integrator.setSave(save);
    integrator.setMinDistance(PARAMS::R / PARAMS::RP );
	int status = integrator.integrate();
    if (status == 2) integrator.backtrack();
	AlgebraicVector<DA> xf = integrator.getXf();
    double final_time = integrator.getT();
    AlgebraicVector<DA> final_state = er3bp_to_secondary(xf, final_time);
    AlgebraicVector<DA> final_spherical = simplified_state_to_oe(final_state, final_time);
    return final_spherical;
}

/* @brief Obtain an element of the Jacobian for a given flow mapping.
 * @param[in] flowMap: Vector of DA objects.
 * @param[in] perturbationVector: Where to locally evaluate the flow map. A vector of zeros implies evaluation at the
 *                                reference.
 * @param[in] i: The row-th element of the Jacobian to compute.
 * @param[in] j: The column-th element of the Jacobian to compute. Indexed from zero!
 * @returns The relevant entry of the Jacobian.
 */
DACE::DA getElementofFlowMap(const DACE::AlgebraicVector<DACE::DA>& flowMap,
                             const DACE::AlgebraicVector<DACE::DA>& perturbationVector,
                             int i, int j){
#ifndef RELEASE
    if (j > DACE::DA::getMaxVariables() || j < 0) throw std::runtime_error("You've requested the derivative with "
                                                                           "respect to a DA variable that does not "
                                                                           "exist.");
    if (i > flowMap.size() || i < 0) throw std::runtime_error("You've requested the derivative with respect to a "
                                                              "non-existent dimension.");
#endif
    DACE::DA ans = flowMap[i].deriv(j+1).eval(perturbationVector);
    return ans;
}

/* @brief Obtain the Cauchy-Green Strain Tensor for a given flow, returning a double-precision representation.
 * @param[in] flowMap: Vector of DA objects.
 * @param[in] deltaK: Where to locally evaluate the flow map. A vector of zeros implies evaluation at the reference.
 * @returns Eigen matrix of double-precision numbers representing the CGST.
 */
Eigen::Matrix<double, 3, 3> getCGST( const DACE::AlgebraicVector<DACE::DA>& flowMap,
                                     const DACE::AlgebraicVector<DACE::DA>& deltaK ) {
#ifndef RELEASE
    if ( flowMap.size() != deltaK.size() || flowMap.empty() || deltaK.empty() ){
        throw std::runtime_error("Invalid inputs to getCGSTDA");
    }
#endif
    namespace e = Eigen;
    e::Matrix<double, 3, 3> cgst;
    for (auto i = 0; i < 3; ++i)
    {
        for (auto j = 0; j < 3; ++j)
        {
            cgst(i,j) = 0.0;
            for (auto k = 0; k < 3; ++k)
            {
                DACE::DA test1 = getElementofFlowMap(flowMap, deltaK, k, i);
                DACE::DA test2 = getElementofFlowMap(flowMap, deltaK, k, j);
                cgst(i, j) += DACE::cons(test1 * test2);
            }
        }
    }
    return cgst;
}

/* @brief Obtains the Cauchy-Green Strain Tensor for a given flow at the reference point.
 * @param[in] flowMap: Vector of DA objects.
 * @returns Eigen matrix of double-precision numbers representing the CGST.
 */
Eigen::Matrix<double, 3, 3> getCGST( const DACE::AlgebraicVector<DACE::DA>& flowMap ) {
#ifndef RELEASE
    if ( flowMap.empty() ){
        throw std::runtime_error("Invalid inputs to getCGSTDA");
    }
#endif
    namespace e = Eigen;
    DACE::AlgebraicVector<DACE::DA> deltaK = {0., 0., 0.};
    e::Matrix<double, 3, 3> cgst;
    for (auto i = 0; i < 3; ++i)
    {
        for (auto j = 0; j < 3; ++j)
        {
            cgst(i,j) = 0.0;
            for (auto k = 0; k < 3; ++k)
            {
                DACE::DA test1 = getElementofFlowMap(flowMap, deltaK, k, i);
                DACE::DA test2 = getElementofFlowMap(flowMap, deltaK, k, j);
                cgst(i, j) += DACE::cons(test1 * test2);
            }
        }
    }
    return cgst;
}

/* @brief Obtain the Cauchy-Green strain tensor for a given flow, returning a DA representation.
 * @param[in] flowMap: Vector of DA objects representing the advected flow.
 * @param[in] deltaK: Where to locally evaluate teh flow map. A vector fo zeros implies evaluation at the reference.
 * @returns Vector-of-vectors representing a matrix of DA objects (i.e. the CGST).
 */
auto getCGSTDA( const std::vector< DACE::DA >& flowMap, const std::vector< DACE::DA >& deltaK )
                                                            -> std::vector< std::vector< DACE::DA > >
{
    using std::vector;
#ifndef RELEASE
	if ( flowMap.size() != deltaK.size() || flowMap.empty() || deltaK.empty() ){
        throw std::runtime_error("Invalid inputs to getCGSTDA");
    }
#endif
	std::vector< std::vector< DACE::DA > > derivative, derivativeTranspose, result;
	/* Give vectors size & zero pad. */
	for (size_t i = 0; i < flowMap.size(); ++i)
	{
		std::vector< DACE::DA > this_row( flowMap.size() );
		for (auto& el : this_row) el = 0.0;
		derivative.emplace_back( this_row );			// Copies - no shared reference between derivative & derivativeTranspose.
		derivativeTranspose.emplace_back( this_row );	// Copies - no shared reference between derivative & derivtiveTranspose.
	}
	/* Give matrices correct entries */
	for (size_t i = 0; i < flowMap.size(); ++i)
	{
		for (size_t j = 0; j < flowMap.size(); ++j)
		{
			derivative[i][j] = flowMap[i].deriv(j+1);
			derivativeTranspose[i][j] = flowMap[j].deriv(i+1);
		}
	}
	/* Compute CGST */
    result = derivativeTranspose * derivative;
	return result;
}

/* @brief Obtain the Cauchy-Green strain tensor for a given flow evaluated at the origin, returning a DA representation.
 * @param[in] flowMap: Vector of DA objects representing the advected flow.
 * @returns Vector-of-vectors representing a matrix of DA objects (i.e. the CGST).
 */
auto getCGSTDA( const std::vector< DACE::DA >& flowMap ) -> std::vector< std::vector< DACE::DA > >
{
   DACE::AlgebraicVector<DACE::DA> deltaK = {0., 0., 0.};
    return getCGSTDA(flowMap, deltaK);
}


/* @brief Obtain the dominant eigenvector of a given Eigen matrix.
 * @param[in] Eigen matrix of double-precision numbers to compute.
 * @returns Eigen vector of double-precision numbers representing the eigenvectors.
 */
template <typename T>
requires std::is_arithmetic_v<T>
auto getDominantEigenvector(const Eigen::Matrix<T, 3, 3>& matrix) -> Eigen::Matrix<double, 3, 1>
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, 3, 3>> solver;
    solver.compute(matrix);
    auto argMax = 0;
    auto eigenvalues = solver.eigenvalues();
    for (auto i = 0; i < eigenvalues.size(); ++i)
    {
        if (eigenvalues[i] > eigenvalues[argMax]) argMax = i;
    }
    Eigen::Matrix<double, 3, 1> eigenvectors = solver.eigenvectors().col(argMax).real();
    return eigenvectors;
}

/* @brief Advect an initial condition forward and obtain its dominant eigenvector.
 * @param[in] x0: Initial condition
 * @param[in] t0: Initial time.
 * @param[in] t: Final time.
 * @returns AlgebraicVector containing the dominant eigenvector.
 */
auto getDominantEigenvector(const DACE::AlgebraicVector<DACE::DA>& x0, const double& t0, double& t)
                            -> DACE::AlgebraicVector<double>
{
	using DACE::DA;
	using Eigen::Matrix;
	using DACE::AlgebraicVector;
	/* Integrate forward, then get CGST, then helicity */
    auto final_condition = getFlowMap(x0, t0, t);
    AlgebraicVector<AlgebraicVector<DA>> jacobian = makeMatrix( 3, 3 );
    for (size_t i = 0; i < final_condition.size(); ++i){
        for (size_t j = 0; j < final_condition.size(); ++j){
            jacobian[i][j] = final_condition[i].deriv(j+1);
        }
    }
    AlgebraicVector<AlgebraicVector<DA>> cgst = transpose(jacobian) * jacobian;
	Matrix<double, 3, 3> cgstDouble;
    for (size_t i = 0; i < 3; ++i){
        for (size_t j = 0; j < 3; ++j){
            cgstDouble(i,j) = DACE::cons(cgst[i][j]);
        }
    }
    Matrix<double, 3, 1> dominantEigenvector = getDominantEigenvector( cgstDouble );
	AlgebraicVector<double> dominantEigenvectorDouble( dominantEigenvector.size() ) ;
	for (auto i = 0; i < dominantEigenvector.size(); ++i) dominantEigenvectorDouble[i] = dominantEigenvector[i];
	return dominantEigenvectorDouble;
}

/* @brief Obtain the helicity for a given initial condition using differential algebra. Note the method is blind to any
 * parameterisation of x0, and acts only on the result of getFlowMap using these initial conditions.
 * @param[in] x0: The initial condition in the given parameterisation.
 * @param[in] t0: The initial time for the integration (i.e. the time at which the initial conditions are valid.
 * @param[in] tf: The final time for the integration (i.e. the time over which the flow map is defined.)
 * @param[inout] helicity: The helicity corresponding to this initial condition and flow time.
 * @returns Eigen::Matrix of doubles containing the 3x3 CGST
 */
auto getHelicityDA(const DACE::AlgebraicVector<DACE::DA>& x0, const double& t0,
                   double& tf, double& helicity) -> Eigen::Matrix<double, 3, 3> {
    using DACE::AlgebraicVector;
    using DACE::DA;
    /* Supplement the initial condition with the DA identity. */
    //auto initial_condition = x0 + AlgebraicVector<DA>::identity();
    /* Flow the initial condition under the mapping. */
    auto final_condition = getFlowMap(x0, t0, tf);
    /* Construct CGST as 6x3 */
    AlgebraicVector<AlgebraicVector<DA>> jacobian = makeMatrix(final_condition.size(), 3);
    for (size_t i = 0; i < final_condition.size(); i++){
        for (size_t j = 0; j < 3; j++){
            jacobian[i][j] = final_condition[i].deriv(j+1);
        }
    }
    /* CGST is jacobian^T * jacobian */
    AlgebraicVector<AlgebraicVector<DA>> cgst = transpose(jacobian) * jacobian;
    /* Construct eigen representation of a double-precision version of the CGST (constant part of the expansion) for use
     * in the eigendecomposition. */
    Eigen::Matrix<double, 3, 3> double_cgst;
    for (size_t i = 0; i < cgst.size(); ++i){
        for (size_t j = 0; j < cgst[i].size(); ++j){
            double_cgst(i, j) = DACE::cons(cgst[i][j]);
        }
    }
    /* Get the double-precision eigenvector. */
    Eigen::Matrix<double, 3, 1> double_ev = getDominantEigenvector(double_cgst);
    AlgebraicVector<DA> dominant_eigenvector( double_ev.size() );
    for (size_t i = 0; i < dominant_eigenvector.size(); ++i){
        dominant_eigenvector[i] = double_ev(i);
    }
    /* Perform power iteration on the double-precision EV to get DA expansion. */
    AlgebraicVector<DA> converged_ev = iterateToConvergence<AlgebraicVector<AlgebraicVector<DA>>,
                                                            AlgebraicVector<DA>>(cgst, dominant_eigenvector);
    /* Construct a Jacobian of the EV */
    Eigen::Matrix<double, 3, 3> ev_jacobian;
    for (size_t i = 0; i < converged_ev.size(); ++i)
    {
        double_ev(i) = DACE::cons(converged_ev[i]);
        for (size_t j = 0; j < 3; ++j){
            ev_jacobian(i, j) = DACE::cons( converged_ev[i].deriv(j+1) );
        }
    }
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(ev_jacobian);
	double cond = svd.singularValues()(0) 
    / svd.singularValues()(svd.singularValues().size()-1);
    /* Now, get the curl. */
    Eigen::Matrix<double, 3, 1> curl;
    curl(0) = ev_jacobian(2,1) - ev_jacobian(1,2); // dfzdy - dfydz
    curl(1) = ev_jacobian(0,2) - ev_jacobian(2,0); // dfxdz - dfzdx
    curl(2) = ev_jacobian(1,0) - ev_jacobian(0,1); // dfydx - dfxdy
    /* Helicity is curl dotted with EV. */
    helicity = fabs( curl.dot(double_ev) );
    return double_cgst;
}                     

#endif //ACROBAT_EIGENVECTORS_H
