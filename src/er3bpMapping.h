/*
 * Created by Jack Tyler on 07/12/2021.
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

#ifndef ACROBAT_ER3BPMAPPING_H
#define ACROBAT_ER3BPMAPPING_H

#include <iostream>
#include "Params.hpp"
#include <vector>
#include <numeric>
#include <algorithm>
#include <dace/dace_s.h>
#include <boost/math/constants/constants.hpp>
#include "linearAlgebra.h"
#include <eigen3/Eigen/Core>
extern "C"
{
#include <SpiceUsr.h>
};

std::ofstream& operator<<(const DACE::AlgebraicVector<DACE::DA>& obj, std::ofstream& out)
{
    for (const auto& el : obj) out << DACE::cons(el) << " ";
    return out;
}

template <typename T>
auto sgn(const T& in) -> int {
	return static_cast<int>(DACE::cons(in) > 0) - static_cast<int>(DACE::cons(in) < 0);
}

/* @brief Construct the Euler rotation matrix about the +z corresponding to a rotation of angle F.
 * @param[in] F: Angle to rotate by.
 * @returns Eigen matrix of the rotation matrix.
 */
auto get_rotation_matrix(const double& F) -> Eigen::Matrix<DACE::DA, 3, 3>
{
    Eigen::Matrix<DACE::DA, 3, 3> r_z;
    /* Eigen will try to construct the matrix of zeros by calling the constructor DA(0), which is equal to 1 in Eigen.
     * Thus, we do it ourselves!
     */
    for (auto i = 0; i < 3; ++i)
    {
        for (auto j = 0; j < 3; ++j)
        {
            r_z(i, j) = 0.0;
        }
    }
    r_z(0,0) = cos(F);
    r_z(0,1) = sin(F);
    r_z(1,0) = -sin(F);
    r_z(1,1) = cos(F);
    r_z(2,2) = 1.0;
    return r_z;
}

/* @brief Construct the derivative of the Euler rotation matrix about the z-axis with an angle F, with respect to F.
 * @param[in] Angle to rotate by.
 * @returns Eigen matrix containing the derivative of R(z).
 */
auto get_rotation_matrix_derivative(const double& F) -> Eigen::Matrix<DACE::DA, 3, 3>
{
    Eigen::Matrix<DACE::DA, 3, 3> r_z_prime;
    /* Eigen will try to construct the matrix of zeros by calling the constructor DA(0), which is equal to 1 in Eigen.
     * Thus, we do it ourselves!
     */
    for (auto i = 0; i < 3; ++i)
    {
        for (auto j = 0; j < 3; ++j)
        {
            r_z_prime(i, j) = 0.0;
        }
    }
    r_z_prime(0,0) = -sin(F);
    r_z_prime(0,1) = cos(F);
    r_z_prime(1,0) = -cos(F);
    r_z_prime(1,1) = -sin(F);
    return r_z_prime;
}

/* @brief Return the time-derivative of the true anomaly F.
 * @param[in] F: True anomaly.
 * @returns Time derivative of the true anomaly.
 */
auto get_f_dot(const double& F)
{
    return ( sqrt(PARAMS::hostGM + PARAMS::targetGM) * (1.0 + PARAMS::EP * cos(F)) * (1.0 + PARAMS::EP * cos(F)) ) /
    ( pow(PARAMS::SMAP, 1.5) * pow((1.0 - PARAMS::EP * PARAMS::EP), 1.5) );
}

/* @brief Return the distance between the target planet and the host planet at a true anomaly F.
 * @param[in] F: True anomaly.
 * @returns The distance between teh target planet and the host planet at the true anomaly.
 */
auto get_sun_mars_distance(const double& F) -> double
{
    return PARAMS::SMAP * (1.0 - PARAMS::EP * PARAMS::EP) / (1.0 + PARAMS::EP * cos(F) );
}

/* @brief Return the F-derivative of the distance between the target planet and the host planet at a true anomaly F.
 * @param[in] F: True anomaly.
 * @returns The derivative of the distance between the host and target planets with respect to F.
 */
auto get_drdf(const double& F) -> double
{
    return -PARAMS::EP * sin(F) / ( PARAMS::SMAP * (1.0 - PARAMS::EP * PARAMS::EP) );
}

/* @brief Obtain the velocity vector associated with a given position.
 * @param[in] state: Cartesian position.
 * @param[in] F: true anomaly
 * @returns Eigen matrix containing the velocity.
 */
auto velocity_from_position(const DACE::AlgebraicVector<DACE::DA>& state, const double& F) -> Eigen::Matrix<DACE::DA, 3, 1>
{
    /* Radius */
    DACE::DA r = sqrt( state[0] * state[0] + state[1] * state[1] + state[2] * state[2] );
    /* Velocity magnitude: comes from OEs */
    DACE::DA velocity_magnitude = sqrt( PARAMS::targetGM * (1 + PARAMS::ECC) / r);
    Eigen::Matrix<DACE::DA, 3, 1> state_eigen;
    for (auto i = 0; i < 3; ++i) state_eigen(i) = state[i];
    Eigen::Matrix<DACE::DA, 3, 1> velocity_vector;
    DACE::AlgebraicVector<DACE::DA> Omega( 3 ), n( 3 ), v_dir(3), velocity_vector_da(3);
    DACE::DA v_dir_norm = 0.0;
    /* Normal vector should be Omega x r */
    Omega[0] = cos(PARAMS::LONGTD); Omega[1] = sin(PARAMS::LONGTD); Omega[2] = 0.0;
    n = sgn(state[2]) * cross3(Omega, state);
	if (DACE::cons(n[2]) < 0.0) n *= -1.;
    DACE::DA n_norm = 0.0, n_norm_sqrt = 0.0;
    for (size_t i = 0; i < n.size(); ++i) n_norm += n[i] * n[i];
    n_norm_sqrt = sqrt(n_norm);
    for (size_t i = 0; i < n.size(); ++i) n[i] /= n_norm_sqrt;
    v_dir = cross3(n, state);
    for (size_t i = 0; i < 3; ++i) v_dir_norm += v_dir[i] * v_dir[i];
    v_dir /= sqrt(v_dir_norm);
    for (size_t i = 0; i < v_dir.size(); ++i) velocity_vector(i) = velocity_magnitude * v_dir[i];
    return velocity_vector;
}

/* @brief Returns the ER3BP velocity from the cartesian state.
 * @param[in] state: Initial cartesian state (inertial frame)
 * @param[in] F: True anomaly at the desired time.
 * @return Eigen matrix containing teh ER3BP velocity from the cartesian state.
 */
auto get_x_prime(const DACE::AlgebraicVector<DACE::DA>& state, const double& F) -> Eigen::Matrix<DACE::DA, 3, 1>
{
    Eigen::Matrix<DACE::DA, 3, 1> dx_dt = velocity_from_position(state, F);
    DACE::DA f_dot = get_f_dot(F);
    return dx_dt / f_dot;
}

/* @brief Transform the ER3BP velocity back into the inertial frame.
 * @param[in] er3bp_velocity: Velocity in the synodic frame (dimensionless)
 * @param[in] inertial_position: Position in the inertial frame (dimensional)
 * @param[in] f: True anomaly.
 * @returns AlgebraicVector containing the inertial velocity.
 */
auto invert_er3bp_velocity(const DACE::AlgebraicVector<DACE::DA>& er3bp_velocity,
                           const DACE::AlgebraicVector<DACE::DA>& inertial_position,
                           double f) -> DACE::AlgebraicVector<DACE::DA>
{
    DACE::DA f_dot = get_f_dot(f);
    DACE::DA radius = get_sun_mars_distance(f);
    Eigen::Matrix<DACE::DA, 3, 3> r_z = get_rotation_matrix(f);
    r_z.transposeInPlace();
    Eigen::Matrix<DACE::DA, 3, 3> r_z_deriv = get_rotation_matrix_derivative(f);
    Eigen::Matrix<DACE::DA, 3, 1> er3bp_velocity_eigen, inertial_position_eigen;
    for (auto i = 0; i < 3; ++i){
        inertial_position_eigen(i) = inertial_position[i];
		er3bp_velocity_eigen(i) = er3bp_velocity[i];
    }
    DACE::DA dr_df = get_drdf(f);
    // See paper for derivation.
    Eigen::Matrix<DACE::DA, 3, 1> v_inertial = f_dot * radius * r_z * er3bp_velocity_eigen -
            f_dot * r_z * r_z_deriv * inertial_position_eigen - f_dot * radius * inertial_position_eigen * dr_df;
    DACE::AlgebraicVector<DACE::DA> final_result(3);
    for (auto i = 0; i < 3; ++i) final_result[i] = v_inertial(i);
    return final_result;
}

template <typename T>
auto simplified_state_to_oe(const DACE::AlgebraicVector<T>& state, double f) -> DACE::AlgebraicVector<T> {
    using DACE::AlgebraicVector;
    AlgebraicVector<T> output_elements;
    AlgebraicVector<T> position( 3 ), velocity( 3 );
    for (size_t i = 0; i < position.size(); ++i){
        position[i] = state[i];
        velocity[i] = state[i+3];
    }
    T r_mag = sqrt( std::inner_product( state.begin(), state.begin() + 3, state.begin(), static_cast<T>(0.0)) );
    AlgebraicVector<T> ang_momentum = cross3(position, velocity);
    T ang_momentum_mag = sqrt( std::inner_product( ang_momentum.begin(), ang_momentum.end(), ang_momentum.begin(),
                                                   static_cast<T>(0.0)) );
    T inc = acos(ang_momentum[2] / ang_momentum_mag);
    T Omega = atan2( ang_momentum[0], -ang_momentum[1] );
    AlgebraicVector<T> normalized_position = position;
    normalizeVector(normalized_position);
    AlgebraicVector<T> eccentricity_vector = cross3(velocity, ang_momentum) / PARAMS::targetGM -
                                             normalized_position;
    T arg_p = atan2(eccentricity_vector[2] / sin(inc), eccentricity_vector[1] * sin(Omega) +
                                                       eccentricity_vector[0] * cos(Omega) );
    output_elements = {r_mag, inc, arg_p};
    return output_elements;
}

template <typename T>
auto spherical_to_cartesian(const DACE::AlgebraicVector<T>& state, double f) -> DACE::AlgebraicVector<DACE::DA> {
    using DACE::DA;
    using DACE::AlgebraicVector;
    DA r = state[0], theta = state[1], phi = state[2];
    DA x = r * cos(theta) * sin(phi);
    DA y = r * sin(theta) * sin(phi);
    DA z = r * cos(phi);
    AlgebraicVector<DA> output = {x, y, z};
    return output;
}

template <typename T>
auto cartesian_to_spherical(const DACE::AlgebraicVector<T>& state, double f) -> DACE::AlgebraicVector<DACE::DA> {
    using DACE::DA;
    using DACE::AlgebraicVector;
    DA x = state[0], y = state[1], z = state[2];
    DA r = sqrt( x * x + y * y + z * z );
    DA th = atan2(y, x);
    DA phi = acos(z / r);
    AlgebraicVector<DA> rthphi = {r, th, phi};
    return rthphi;
}


/* @brief Convert a state (position, velocity) from the ER3BP synodic frame and units back to the inertial frame.
 * @param[in] er3bp: ER3BP state
 * @param[in] F: true anomaly at the desired time.
 * @returns AlgebraicVector containing the inertial state.
 */
auto er3bp_to_secondary(const DACE::AlgebraicVector<DACE::DA>& er3bp, double f) -> DACE::AlgebraicVector<DACE::DA>
{
    using DACE::AlgebraicVector;
    const double er3bp_mu = PARAMS::targetGM / (PARAMS::hostGM + PARAMS::targetGM);
    Eigen::Matrix<DACE::DA, 3, 3> r_z = get_rotation_matrix(f);
    Eigen::Matrix<DACE::DA, 3, 1> er3bp_eigen, translation_vector, inertial_eigen;
    r_z.transposeInPlace();
    for (auto i = 0; i < 3; ++i) er3bp_eigen(i) = er3bp[i];
    translation_vector(0) = 1.0; translation_vector(1) = 0.0; translation_vector(2) = 0.0;
    DACE::DA radius = get_sun_mars_distance(f);
    inertial_eigen = r_z * (er3bp_eigen - (1.0 - er3bp_mu)  * translation_vector);
    inertial_eigen *= radius;
    DACE::AlgebraicVector<DACE::DA> inertial_position(3), er3bp_velocity(3);
    for (auto i = 0; i < 3; ++i){
        inertial_position[i] = inertial_eigen(i);
        er3bp_velocity[i] = er3bp[i+3];
    }
    DACE::AlgebraicVector<DACE::DA> inertial_velocity = invert_er3bp_velocity(er3bp_velocity, inertial_position, f);
    DACE::AlgebraicVector<DACE::DA> inertial_state;
    std::copy( inertial_position.begin(), inertial_position.end(), std::back_inserter(inertial_state) );
    std::copy( inertial_velocity.begin(), inertial_velocity.end(), std::back_inserter(inertial_state) );
    return inertial_state;
}

void oes_to_state_simplified(const DACE::AlgebraicVector<DACE::DA>& oes, double& mu, DACE::AlgebraicVector<DACE::DA>& state)
{
    using DACE::DA;
    DA a = oes[0] / (1. - oes[1]);
    DA e = oes[1];
    DA i = oes[2];
    DA omg = oes[3];
    DA omp = oes[4];
    DA M = oes[5];
    DA R[3][3];
    /* If hyperbolic, switch to negative sma */
    if (DACE::cons(e) > 1) a *= -1;
    /* Get position, velocity in the perifocal frame */
    DA xper, yper, xdotper, ydotper;
    /* Switch from mean anomaly to eccentric anomaly */
    DA EA = M; // Always start at zero mean anomaly anyway!
    if (DACE::cons(e) < 1.)
    {
        DA b, n;
        b = a * sqrt(1 - e * e);
        n = sqrt(mu) / sqrt( a * a *  a);

        xper = a * (cos(EA) - e);
        yper = b * sin(EA);
        xdotper = -(a * n * sin(EA)) / (1 - e * cos(EA));
        ydotper = (b * n * cos(EA)) / (1 - e * cos(EA));
    } else // EA is the Gudermannian
    {
        DA b, n;
        b = -a * sqrt(e * e - 1);
        n = sqrt(-mu / (a * a * a));

        DA dNdZeta = e * (1 + tan(EA) * tan(EA)) - (0.5 + 0.5 * pow(tan(0.5 * EA + M_PI_4), 2)) / tan(0.5 * EA + M_PI_4);

        xper = a / cos(EA) - a * e;
        yper = b * tan(EA);

        xdotper = a * tan(EA) / cos(EA) * n / dNdZeta;
        ydotper = b / pow(cos(EA), 2) * n / dNdZeta;
    }
    // 2 - We then built the rotation matrix from perifocal reference frame to inertial

    DA cosomg = cos(omg);
    DA cosomp = cos(omp);
    DA sinomg = sin(omg);
    DA sinomp = sin(omp);
    DA cosi = cos(i);
    DA sini = sin(i);

    R[0][0] = cosomg * cosomp - sinomg * sinomp * cosi;
    R[0][1] = -cosomg * sinomp - sinomg * cosomp * cosi;
    R[0][2] = sinomg * sini;
    R[1][0] = sinomg * cosomp + cosomg * sinomp * cosi;
    R[1][1] = -sinomg * sinomp + cosomg * cosomp * cosi;
    R[1][2] = -cosomg * sini;
    R[2][0] = sinomp * sini;
    R[2][1] = cosomp * sini;
    R[2][2] = cosi;
    DA temp[3] = {xper, yper, 0.0};
    DA temp2[3] = {xdotper, ydotper, 0};
    DA r0[3], v0[3];
    for (int j = 0; j < 3; j++) {
        r0[j] = 0.0;
        v0[j] = 0.0;
        for (int k = 0; k < 3; k++) {
            r0[j] += R[j][k] * temp[k];
            v0[j] += R[j][k] * temp2[k];
        }
    }
    for (int j = 0; j < 3; ++j)
    {
        state[j] = r0[j];
        state[j+3] = v0[j];
    }
}

/* @brief Transform a position from the inertial frame to the ER3BP.
 * @param[in] Inertial: state in dimensional, inertial units
 * @param[in] f: True anomaly at the time associated with inertial.
 * @returns AlgebraicVector containing the ER3BP state.
 */
auto manual_position_vector_er3bp(const DACE::AlgebraicVector<DACE::DA>& inertial,
                                  const double& f) -> DACE::AlgebraicVector<DACE::DA>
{
	using DACE::DA;
	Eigen::Matrix<DA, 3, 3> r_z = get_rotation_matrix(f);
	Eigen::Matrix<DA, 3, 1> inertial_eigen;
	for (auto i = 0; i < 3; ++i) inertial_eigen[i] = inertial[i];
	DA distance = get_sun_mars_distance(f);
	Eigen::Matrix<DA, 3, 1> translation_vector;	
	translation_vector(0) = 1.0 - PARAMS::MU; translation_vector(1) = 0.0; translation_vector(2) = 0;
	Eigen::Matrix<DA, 3, 1> r_rp_eigen;
	r_rp_eigen = r_z * inertial_eigen / distance + translation_vector;
	DACE::AlgebraicVector<DA> output(3);
	for (auto i = 0; i < 3; ++i) output[i] = r_rp_eigen[i];
	return output;
}

/* @brief Transform a position from the inertial frame to the ER3BP.
 * @param[in] State: position in the cartesian inertial frame.
 * @param[in] velocity: velocity in the cartesian frame.
 * @param[in] F: true anomaly at the desired time.
 * @returns The velocity in the ER3BP frame.
 */
auto manual_velocity_vector_test(const DACE::AlgebraicVector<DACE::DA>& state,
                                 const DACE::AlgebraicVector<DACE::DA>& velocity,
                                 const double& F) -> DACE::AlgebraicVector<DACE::DA>
{
    using DACE::DA;
	DA radius = get_sun_mars_distance(F);
	DA f_dot = get_f_dot(F);
	DA drdf = get_drdf(F);
	Eigen::Matrix<DA, 3, 3> dr_zdf = get_rotation_matrix_derivative(F);
	Eigen::Matrix<DA, 3, 3> r_z = get_rotation_matrix(F);
	Eigen::Matrix<DA, 3, 1> v_p, state_eigen, velocity_eigen, rotation_derivative0;
	for (auto i = 0; i < 3; ++i){
        state_eigen(i) = state[i];
        velocity_eigen(i) = velocity[i];
    }
    v_p = dr_zdf * state_eigen / radius + r_z * state_eigen * drdf + r_z * velocity_eigen / (radius * f_dot);
	DACE::AlgebraicVector<DA> result(3);
	for (auto i = 0; i < 3; ++i) result[i] = v_p(i);
	return result;
}

template <typename T> requires arithmetic<T>
auto simplified_oe_to_state(const DACE::AlgebraicVector<T>& secondary, double f) -> DACE::AlgebraicVector<T> {
    /* See paper for derivation - secondary should be the position in parameter space (r, inc, om) */
    using DACE::DA;
    using DACE::AlgebraicVector;
    T radius = secondary[0];
    T inc    = secondary[1];
    T arg_p  = secondary[2];
    T velocity = sqrt(PARAMS::targetGM * (1 + PARAMS::ECC) / radius);
    /* Pre-compute relevant terms */
    T sin_argp = sin(arg_p);
    T cos_argp = cos(arg_p);
    T sin_inc  = sin(inc);
    T cos_inc  = cos(inc);
    AlgebraicVector<T> state( 6 );     // Hardcode for 6D phase space
    state[0] = radius * cos_argp;
    state[1] = radius * sin_argp * cos_inc;
    state[2] = radius * sin_argp * sin_inc;
    state[3] = velocity * -sin_argp;
    state[4] = velocity * cos_argp * cos_inc;
    state[5] = velocity * cos_argp * sin_inc;
    return state;
}

/* @brief Transform a state (position, velocity) from the inertial frame to the ER3BP frame.
 * @param[in] secondary: State in the target-centered inertial frame.
 * @param[in] f: True anomaly.
 * @returns AlgebraicVector containing the ER3BP state.
 */
auto secondary_to_er3bp(const DACE::AlgebraicVector<DACE::DA>& secondary, double f) -> DACE::AlgebraicVector<DACE::DA>
{
    using DACE::AlgebraicVector;
    using DACE::DA;
	AlgebraicVector<DA> eme_state = secondary;
    AlgebraicVector<DA> eme_position = {eme_state[0], eme_state[1], eme_state[2]};
    AlgebraicVector<DA> eme_velocity = {eme_state[3], eme_state[4], eme_state[5]};
    AlgebraicVector<DA> er3bp_position = manual_position_vector_er3bp(eme_position, f);
    AlgebraicVector<DA> er3bp_velocity = manual_velocity_vector_test(eme_position, eme_velocity, f);
    AlgebraicVector<DA> er3bp_state;
    std::copy( er3bp_position.begin(), er3bp_position.end(), std::back_inserter(er3bp_state) );
    std::copy( er3bp_velocity.begin(), er3bp_velocity.end(), std::back_inserter(er3bp_state) );
    return er3bp_state;
}

#endif //ACROBAT_ER3BPMAPPING_H
