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

#ifndef ACROBAT_INTEGRATION_H
#define ACROBAT_INTEGRATION_H

#include <vector>
#include <boost/numeric/odeint.hpp>
#include <dace/dace_s.h>
#include "rhs.h"
#include <fstream>
#include "Params.hpp"

/* Adds interoperability between ODEINT and DACE */
namespace boost{ namespace numeric{ namespace odeint{
            /* @brief Compute absolute value of an AlgebraicVector (abs of each element separately.)
             * @param[in] p: Vector to abs.
             * @returns AlgebraicVector containing the abs of p.
             */
            template <typename T>
            DACE::AlgebraicVector<T> abs(DACE::AlgebraicVector<T>& p)
            {
                DACE::AlgebraicVector<T> res( p.size() );
                for (auto i = 0; i < res.size(); ++i)
                {
                    res[i] = abs( p[i] );
                }
                return res;
            }
            /* @brief Compute absolute value of a const AlgebraicVector (abs of each element separately.)
             * @param[in] p: Vector to abs.
             * @returns AlgebraicVector containing the abs of p.
             */
            template <typename T>
            DACE::AlgebraicVector<T> abs(const DACE::AlgebraicVector<T>& p)
            {
                DACE::AlgebraicVector<T> res( p.size() );
                for (auto i = 0; i < res.size(); ++i)
                {
                    res[i] = abs( p[i] );
                }
                return res;
            }
            /* @brief Template specialisation for the L2 norm of an AlgebraicVector */
            template <>
            struct vector_space_norm_inf<DACE::AlgebraicVector<DACE::DA>>
            {
                typedef double result_type;
                double operator()(const DACE::AlgebraicVector<DACE::DA>& p) const
                {
                    double maxCoeff = 0.0;
                    for (unsigned i = 0; i < p.size(); ++i) maxCoeff = std::max( abs( p[i] ) , maxCoeff );
                    return maxCoeff;
                }
            };
            template <>
            struct vector_space_norm_inf<DACE::AlgebraicVector<double>>
            {
                typedef double result_type;
                double operator()(const DACE::AlgebraicVector<double>& p) const
                {
                    double maxCoeff = 0.0;
                    for (unsigned i = 0; i < p.size(); ++i) maxCoeff = std::max( std::abs( DACE::cons(p[i]) ), maxCoeff );
                    return maxCoeff;
                }
            };
            template <>
            struct vector_space_norm_inf<std::vector<DACE::DA>>
            {
                typedef double result_type;
                double operator()(const std::vector<DACE::DA>& p) const
                {
                    double maxCoeff = 0.0;
                    for (size_t i = 0; i < p.size(); ++i) maxCoeff = std::max( abs( p[i]) , maxCoeff );
                    return maxCoeff;
                }
            };

            template<> struct is_resizeable<DACE::AlgebraicVector<DACE::DA>>
            {
                typedef boost::false_type type;
                const static bool value = type::value;
            };
            template<> struct is_resizeable<DACE::AlgebraicVector<double>>
            {
                typedef boost::false_type type;
                const static bool value = type::value;
            };
        }
    }
}

template <class vectorType, class timeType>
class integration
{
public:
    integration() = default;
    integration( DACE::AlgebraicVector<vectorType> x0, timeType t ) : _x0(x0), _xf(x0), _t0(0.0), _tf(t), _t(_t0) {};
    integration( DACE::AlgebraicVector<vectorType> x0, timeType t0, timeType tf ) : _x0(x0), _xf(x0), _t0(t0), _tf(tf),
                                                                                    _t(t0){};

    int integrate() {
	const double er3bp_mu = PARAMS::MU;
        _integrationStatus = 0;
        /* Take a step, compute the position/time at this step and then return this status. */
        auto stepper = boost::numeric::odeint::make_controlled(_tol, _tol, _rk78);
		std::ofstream output;
		output.precision(16);
		if (this->_save){
			std::string this_fname = "trajectory_" + std::to_string(this->_id);
			output.open(this_fname);
		}
        while (_integrationStatus == 0)
        {
            boost::numeric::odeint::controlled_step_result stepSuccess = boost::numeric::odeint::fail;
            while (stepSuccess == boost::numeric::odeint::fail)
            { // Keep attempting to step until success
                stepSuccess = stepper.try_step( this->forceFunction, _xf, _t, _dt);
            }
            numSteps++;
			if (this->_save){
				for (int i = 0; i < this->_xf.size() - 1; ++i){
					output << DACE::cons( this->_xf[i] ) << ",";
				}
				output << DACE::cons( this->_xf[ this->_xf.size() - 1 ] ) << "\n";
			}
            /* Now check the integration statuses */
            timeType distance = sqrt( DACE::cons((_xf[0]-1.+er3bp_mu) * (_xf[0]-1.0+er3bp_mu) + _xf[1]*_xf[1] +
                    _xf[2]*_xf[2]) );
            if (distance <= _minDistance) _integrationStatus = 1;
            if (_t >= _tf) _integrationStatus = 2;
        }
		if (this->_save){
			output.close();
		}
        return _integrationStatus;
    }

    /* @brief From a current time >= tf, go backwards to time tf. */
    void backtrack()
    {
        int integrationDirection = static_cast<int>( (_t - _tf) > 0 ) - static_cast<int>( (_t - _tf) < 0 ); // Sign of integration to perform
        timeType local_dt = _dt * integrationDirection;             // Flip the current direction of integration
        boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(_tol, _tol, _rk78),
                                                   this->forceFunction, _xf, _t, _tf, local_dt);
        _t = _tf;
    }

    /* Getters and setters */
    const DACE::AlgebraicVector <vectorType> &getX0() const {
        return _x0;
    }

    const DACE::AlgebraicVector<vectorType> &getXf() const {
        return _xf;
    }

    void setXf(const DACE::AlgebraicVector<vectorType> &xf) {
        _xf = xf;
    }

    timeType getT0() const {
        return _t0;
    }

    void setT0(timeType t0) {
        _t0 = t0;
    }

    timeType getTf() const {
        return _tf;
    }

    void setTf(timeType tf) {
        _tf = tf;
    }

    const timeType getTol() const {
        return _tol;
    }

    void setTol(timeType tol)
    {
        _tol = tol;
    }

    timeType getMinDistance() const {
        return _minDistance;
    }

    void setMinDistance(timeType minDistance){
        _minDistance = minDistance;
    }

    timeType getT() const {
        return _t;
    }

    void setT(timeType t) {
        _t = t;
    }

    unsigned long getNumSteps() const {
        return numSteps;
    }

    void setNumSteps(unsigned long numSteps) {
        integration::numSteps = numSteps;
    }

    void setForceFunction(
            void (*forceFunction)(const DACE::AlgebraicVector<vectorType> &, DACE::AlgebraicVector<vectorType> &, timeType)) {
        integration::forceFunction = forceFunction;
    }

    void setSave(bool save)
    {
        this->_save = save;
    }

	void setId(int id){
		this->_id = id;
	}

	int getId() const {
		return this->_id;
	}

private:
    DACE::AlgebraicVector<vectorType> _x0;
    DACE::AlgebraicVector<vectorType> _xf;
    timeType _t0;
    timeType _tf;
    timeType _t;
    timeType _dt = 0.1;
    uint_fast8_t _integrationStatus = 0;
    timeType _tol = 1e-13;
    timeType _minDistance = 1e-010;
    unsigned long numSteps = 0;
	int _id = 0;
    bool _save = false;

    boost::numeric::odeint::runge_kutta_fehlberg78< DACE::AlgebraicVector<vectorType>, timeType, DACE::AlgebraicVector<vectorType>,
            timeType, boost::numeric::odeint::vector_space_algebra > _rk78;
    void (*forceFunction)(const DACE::AlgebraicVector<vectorType>&, DACE::AlgebraicVector<vectorType>&, timeType)
    = RHS::forceFunction;
    void (*observer)(const DACE::AlgebraicVector<vectorType>&, timeType) = RHS::forceFunctionObserver;
};

#endif //ACROBAT_INTEGRATION_H
