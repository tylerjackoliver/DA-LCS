/*
 * Created by Jack Tyler on 06/07/2021.
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

#ifndef ACROBAT_HELICITYINTEGRATOR_H
#define ACROBAT_HELICITYINTEGRATOR_H

#include "rhs.h"
#include "integration.h"
#include <vector>
#include <algorithm>
#include <dace/dace.h>
#include "eigenvectors.h"
#include "powerLaw.h"
#include "Params.hpp"
#include <eigen3/Eigen/Core>
#include <boost/numeric/odeint.hpp>

/* @brief Class that performs an integration of the strainlines for the LCS. */
template <class vector, class numeric>
        requires has_size<vector>
class helicityIntegrator
{
public:
    /* @brief Default constructor. */
    helicityIntegrator() = default;
    /* @brief Constructor for a helicityIntegrator class.
     * @param[in] x0: Initial position of the line (seedpoint)
     * @param[in] Threshold: point at which to stop integrating the strainline.
     */
    helicityIntegrator( vector x0, numeric threshold ) : _x0(x0), _x(x0), _threshold(threshold) {};

    /* @brief Force function for the strainline integration.
     * @param[in] x0: Position at this step (either cartesian or parameter space)
     * @param[in] x: Derivative of the state at this step (either cartesian or parameter space)
     * @param[in] t: Current time associated with this step.
     */
    void strainlineRHS( const vector& x0, vector& x, numeric t )
    {
        using DACE::DA;
        using DACE::AlgebraicVector;
        if (x.size() != x0.size()) x.resize( x0.size() );
        /* Integrate state forward */
        vector x0Copy = x0;
        double integrationTime = PARAMS::T0 + PARAMS::Tf;
        auto final_condition = getFlowMap(x0Copy, PARAMS::T0, integrationTime);
        /* Construct CGST in DA and double precision to get dominant eigenvector */
        AlgebraicVector<AlgebraicVector<DA>> jacobian = makeMatrix( 3, 3 );
        for (size_t i = 0; i < 3; ++i){
            for (size_t j = 0; j < 3; ++j){
                jacobian[i][j] = final_condition[i].deriv(j+1);
            }
        }
        AlgebraicVector<AlgebraicVector<DA>> cgst_da = transpose(jacobian) * jacobian;
        Eigen::Matrix<double, 3, 3> cgstDouble;
        for (auto i = 0; i < 3; ++i)
        {
            for (auto j = 0; j < 3; ++j) cgstDouble(i, j) = DACE::cons( cgst_da[i][j] );
        }
        /* Get eigenvector */
        Eigen::Matrix<double, 3, 1> dominantEigenvectorDouble = getDominantEigenvector(cgstDouble);
        /* Define derivative */
        AlgebraicVector<DACE::DA> unitNormal = {0., 1., 0.}; // Normal to +inc-direction
        AlgebraicVector<DACE::DA> dominantEigenvectorVector = {dominantEigenvectorDouble[0],
                                                               dominantEigenvectorDouble[1],
                                                               dominantEigenvectorDouble[2]};
        x = unitNormal.cross( dominantEigenvectorVector );
        double innerProduct = DACE::cons( x.dot( this->_previousSolution ) );
        double dir = innerProduct > 0. ? 1. : -1.;
        for (size_t i = 0; i < x.size(); ++i) x[i] *= dir;
    }

    /* @brief Integrate a single strainline. Writes the output to a file. */
    auto integrateStrainline() -> void
    {
        /* The eigenvector is only defined up to a sign. Thus, go in both directions, just to be sure! */
        for (int dir = -1; dir < 2; dir += 2) {
            /* Instantiate stepper */
            boost::numeric::odeint::runge_kutta_dopri5<vector, numeric, vector, numeric,
                    boost::numeric::odeint::vector_space_algebra> rk45;
            auto stepper = boost::numeric::odeint::make_controlled(this->_integrationTolerance,
                                                                   this->_integrationTolerance,
                                                                   rk45);
            /* Integration variables */
            int integrationStatus = 0;
            double dt = 1e-010;
            DACE::AlgebraicVector<DACE::DA> tmpCondition = this->_x0;
            this->_x = this->_x0;
            double integrationTime = PARAMS::T0 + PARAMS::Tf;
            DACE::AlgebraicVector<double> dominantEigenvector = getDominantEigenvector(tmpCondition, PARAMS::T0,
                                                                                       integrationTime);
            this->_previousSolution = {dir * -dominantEigenvector[1], 0, dir * dominantEigenvector[2]};
            /* Get the helicity at the initial point */
            getHelicityDA(tmpCondition, PARAMS::T0, integrationTime, this->_helicity);
            /* Record the number of steps we've taken and the current running average of helicity. */
            this->_numSteps = 1;
            this->_helicityAvg = this->_helicity / this->_numSteps;
            this->_helicityAverages.template emplace_back(this->_helicityAvg);
            this->_locationHistory.template emplace_back(this->_x);
            /* Now integrate */
            while (integrationStatus == 0) {
                boost::numeric::odeint::controlled_step_result stepResult = boost::numeric::odeint::fail;
                while (stepResult == boost::numeric::odeint::fail) { // Keep trying to step until timestep selected is OK
                    stepResult = stepper.try_step([this](auto const &x, auto &dxdt, auto t) {
                        this->strainlineRHS(x, dxdt, t);
                    }, this->_x, this->_t, dt); // Requires Lambda since it's a member function of a class.
                }
                this->_numSteps++;
                /* Now, we need the helicity at this new point. */
                double thisHelicity = std::numeric_limits<double>::max();
                tmpCondition = this->_x;
                integrationTime = PARAMS::T0 + PARAMS::Tf;
                Eigen::Matrix<double, 3, 3> cgst = getHelicityDA(tmpCondition, PARAMS::T0, integrationTime,
                                                                  thisHelicity);
                /* Now, we need the solution at this point in order to ensure 'smooth' vector field. */
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3> > es;
                es.compute(cgst);
                int argMax = 0;
                /* Now we know all the eigenvectors, we want just the maximum one; I found no way that Eigen
                 * automatically does this, so let's try them all! */
                for (size_t i = 0; i < es.eigenvalues().size(); ++i) {
                    if (es.eigenvalues()[i] > es.eigenvalues()[argMax]) argMax = i;
                }
                DACE::AlgebraicVector<DACE::DA> dominantEigenvectorDA(3);
                for (size_t i = 0; i < dominantEigenvector.size(); ++i) {
                    dominantEigenvectorDA[i] = es.eigenvectors().col(argMax)[i];
                }
                /* Now unit normal, get solution */
                DACE::AlgebraicVector<DACE::DA> unitNormal = {0., 1., 0.};
                DACE::AlgebraicVector<DACE::DA> newSolution = unitNormal.cross(dominantEigenvectorDA);
                double newSign = DACE::cons(newSolution.dot(this->_previousSolution)) > 0 ? 1. : -1.;
                for (auto i = 0; i < 3; ++i) {
                    newSolution[i] *= newSign;
                }
                /* Update results arrays. */
                this->_previousSolution = newSolution;
                this->_helicity += thisHelicity;
                this->_helicityAvg = this->_helicity / this->_numSteps;
                this->_helicityAverages.template emplace_back(this->_helicityAvg);
                this->_locationHistory.template emplace_back(this->_x);
                /* Check stopping conditions */
                if (this->_helicityAvg > this->_threshold) integrationStatus = 1;
                if (this->_numSteps > this->_maxSteps) integrationStatus = 2;
                if ( DACE::cons(this->_x[0]) > PARAMS::RS || DACE::cons(this->_x[0]) < PARAMS::R) integrationStatus = 3;
            }
            /* If this is the first direction, we need to:
             * - reset the running helicity; reset the number of steps;
             * - reset the integration status;
             * - 'reverse' the trajectory to get one smooth line.
             * We also remove the 'last' point in the reversed trajectory, since this is the 'start' point that will be
             * written by the next one.
             */
            if (dir < 0) {
                this->_helicity = 0;
                this->_helicityAvg = 0;
                this->_numSteps = 0;
                std::reverse(this->_locationHistory.begin(), this->_locationHistory.end());
                std::reverse(this->_helicityAverages.begin(), this->_helicityAverages.end());
                this->_locationHistory.pop_back();
                this->_helicityAverages.pop_back();
            }
        }
    }

    /* @brief Write a single strainline out to a file.
     * @param[in] id: This strainline's integer ID, which controls the file naming.
     * File structure is <strainline length> x 4: first three columns are the position (cartesian/parameter space);
     * the remaining column is the running average of helicity.
     */
    void writeStrainline( int id )
    {
        std::ofstream output;
        output.open("strainline_" + std::to_string(id) );
        for (auto i = 0; i < this->_locationHistory.size(); ++i)
        {
            for (auto j = 0; j < this->_locationHistory[i].size(); ++j)
            {
                output << DACE::cons(this->_locationHistory[i][j]) << " ";
            }
            output << DACE::cons( this->_helicityAverages[i] ) << "\n";
        }
        output.close();
    }

    /* @brief Write a single strainline out to a file.
     * @param[in] id: This strainline's integer ID, which controls the file naming.
     * @param[in] path: Prefix to the folder to write the file in.
     * File structure is <strainline length> x 4: first three columns are the position (cartesian/parameter space);
     * the remaining column is the running average of helicity.
     */
    void writeStrainline( int id, const std::string& path )
    {
        std::ofstream output;
        output.open(path + "/strainline_" + std::to_string(id) );
        for (size_t i = 0; i < this->_locationHistory.size(); ++i)
        {
            for (size_t j = 0; j < this->_locationHistory[i].size(); ++j)
            {
                output << DACE::cons(this->_locationHistory[i][j]) << " ";
            }
            output << DACE::cons( this->_helicityAverages[i] ) << "\n";
        }
        output.close();
    }

    /* Getters and setters */
    vector getX0() const {
        return _x0;
    }

    void setX0(vector x0) {
        _x0 = x0;
    }

    vector getX() const {
        return _x;
    }

    void setX(vector x) {
        _x = x;
    }

    numeric getThreshold() const {
        return _threshold;
    }

    void setThreshold(numeric threshold) {
        _threshold = threshold;
    }

    numeric getT() const {
        return _t;
    }

    void setT(numeric t) {
        _t = t;
    }

    numeric getIntegrationTolerance() const {
        return _integrationTolerance;
    }

    void setIntegrationTolerance(numeric integrationTolerance) {
        _integrationTolerance = integrationTolerance;
    }

    numeric getHelicity() const {
        return _helicity;
    }

    void setHelicity(numeric helicity) {
        _helicity = helicity;
    }

    const std::vector<vector> &getLocationHistory() const {
        return _locationHistory;
    }

    void setLocationHistory(const std::vector<vector> &locationHistory) {
        _locationHistory = locationHistory;
    }

    numeric getHelicityAvg() const {
        return _helicityAvg;
    }

    void setHelicityAvg(numeric helicityAvg) {
        _helicityAvg = helicityAvg;
    }

    int getNumSteps() const {
        return _numSteps;
    }

    void setNumSteps(int numSteps) {
        _numSteps = numSteps;
    }

    vector getPreviousSolution() const {
        return _previousSolution;
    }

    void setPreviousSolution(vector previousSolution) {
        _previousSolution = previousSolution;
    }

    numeric getMinX() const {
        return minX;
    }

    void setMinX(numeric minX) {
        helicityIntegrator::minX = minX;
    }

    numeric getMaxX() const {
        return maxX;
    }

    void setMaxX(numeric maxX) {
        helicityIntegrator::maxX = maxX;
    }

    numeric getMinY() const {
        return minY;
    }

    void setMinY(numeric minY) {
        helicityIntegrator::minY = minY;
    }

    numeric getMaxY() const {
        return maxY;
    }

    void setMaxY(numeric maxY) {
        helicityIntegrator::maxY = maxY;
    }

    numeric getMinZ() const {
        return minZ;
    }

    void setMinZ(numeric minZ) {
        helicityIntegrator::minZ = minZ;
    }

    numeric getMaxZ() const {
        return maxZ;
    }

    void setMaxZ(numeric maxZ) {
        helicityIntegrator::maxZ = maxZ;
    }

    int getMaxSteps() const {
        return _maxSteps;
    }

    void setMaxSteps(int maxSteps) {
        _maxSteps = maxSteps;
    }

    void setFudgeFactor(numeric fudgeFactor){
        this->fudgeFactor = fudgeFactor;
    }

private:
    std::vector<vector> _locationHistory;
    std::vector<numeric> _helicityAverages;
    vector _x0, _x, _previousSolution;
    numeric _threshold, _t, _integrationTolerance, _helicity = 0, _helicityAvg = 0;
    numeric minX = PARAMS::R, maxX = PARAMS::RS, minY = 0.0, maxY = 8.0 * atan(1.0), minZ = 0.0, maxZ = 8.0 * atan(1.0);
    numeric fudgeFactor = 2. * 4.0 * atan(1.0);
    int _numSteps = 0, _maxSteps = 1000;

};

#endif //ACROBAT_HELICITYINTEGRATOR_H
