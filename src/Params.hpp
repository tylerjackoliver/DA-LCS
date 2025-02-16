/* 
 * This file is part of the TRANSPORT-BARRIERS distribution (https://github.com/tylerjackoliver/transport-barriers/).
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __ACROBAT_PARAMS_H__
#define __ACROBAT_PARAMS_H__

#include <unordered_map>
#include <string>
#include <vector>
#include <cmath>
#include <numbers>
extern "C"
{
#include <SpiceUsr.h>
}

/* @brief Namespace that contains Parameters used in the simulation. */
namespace PARAMS {
    const int xGridSize = 500;
    const int yGridSize = 500;
    const int zGridSize = 1;

    double integrationAnomaly = 10.0;

    /* Orbital parameters */
    const double ECC = 0.95;
    const double LONGTD = 118 * std::numbers::pi_v<double> / 180.;
    double EPOCH = 649684869.1832148;
    double T0 = 0.0;
    double Tf = 2.0 * std::numbers::pi_v<double>;
    double INC = 35 * std::numbers::pi_v<double> / 180.;
    double EP = 0.0935;
    double M1 = 0.0;
    double M2 = 0.0;
    double MU;
    double SMAP = 2.279e8;
    double helicityThreshold = 1e-018;
    double helicityIntegrationThreshold = 1e-018;

    const double G = 6.6743015e-20; // Gravitational constant km3/s^2

    /* Bodies */
    const std::string TARGET = "Mars";
    const std::string HOST = "Sun";

    const double PI = 4.0 * atan(1.0);
    const double AU = 1.495978707e8;

    /* ~~~~~~ Derived parameters ~~~~~~ */
    /* GM of TARGET */
    double targetGM = 42828.372; // E.g. Earth
    /* GM of the HOST */
    double hostGM = 1.32712440042e11; // E.g. Sun
    /* Planetary radius of the TARGET */
    double R = 3389.5; // Planetary radius, in km
    /* Sphere of influence of the TARGET */
    double RS = 312.81 * R; // Sphere of influence
    /* Maximum integration time for determining a trajectory to be acrobatic */
    double RP = 0.0;
	double RA = 0.0;
    /* GMs of other planets */
    double MUP = 0.0;

	double _RP = RP;
	double _RA = RA;
	double _RS = RS;
	double _R = R;
	double _MUP;
	double _MU;
	double _SMAP;
	double _hostGM;
	double _targetGM;
	
	DACE::AlgebraicVector<DACE::DA> plane_normal;

    std::unordered_map<std::string, double> PLANETARY_GM    ({{"Sun", 1.327e11},
                                                              {"Mercury", 2.203e4},
                                                              {"Jupiter", 1.267e8},
                                                              {"Saturn", 3.794e7},
                                                              {"Pluto", 8.696e2},
                                                              {"Charon", 1.058e2},
                                                              {"Io", 5.960e3},
                                                              {"Callisto", 7.179e3},
                                                              {"Ganymede", 9.888e3},
                                                              {"Triton", 1.428E+03},
                                                              {"Titan", 8.978E+03},
                                                              {"Saturn", 3.794E+07},
                                                              {"Europa", 3.203E+03},
                                                              {"Venus", 3.249e5},
                                                              {"Moon", 4.903e3},
                                                              {"Earth",  3.986e5},
                                                              {"Mars", 42828.372}});
	std::unordered_map<std::string, double> _PLANETARY_GM = PLANETARY_GM;

    /* @brief Initialise run-time parameters: ER3BP mass parameter, planetary gravitational parameters and periapsis
     * radius. Results are stored in the PARAMS:: namespace.
     */
    void initialiseVariables() {
        furnsh_c("../spice/metakernel.tm");
        PARAMS::MUP = PARAMS::hostGM;
        PARAMS::MU = PARAMS::targetGM / (PARAMS::targetGM + PARAMS::hostGM); // 3BP mu
        PARAMS::M1 = PARAMS::hostGM   / PARAMS::G;
        PARAMS::M2 = PARAMS::targetGM / PARAMS::G;
        PARAMS::RP = PARAMS::SMAP * (1.0 - PARAMS::EP);
		PARAMS::RA = PARAMS::SMAP * (1.0 + PARAMS::EP);
		_MUP = MUP;
		_R = R;
		_RS = RS;
		_SMAP = SMAP;
		_MU = MU;
		_RP = RP;
		_RA = RA;
		_PLANETARY_GM = PLANETARY_GM;
		_hostGM = hostGM;
		_targetGM = targetGM;
    }

    /* @brief Normalize all the variables used in the runtime. */
    void normalizeVariables() {
        auto normalizing_radius = PARAMS::RP;
        // Normalize the distance
        PARAMS::RP /= normalizing_radius;
        PARAMS::SMAP /= normalizing_radius;
        PARAMS::R /= normalizing_radius;
        PARAMS::RS /= normalizing_radius;
        // Gravitational parameter
        PARAMS::MUP /= pow( normalizing_radius, 3. );
        PARAMS::targetGM /= pow( normalizing_radius, 3. );
        PARAMS::hostGM /= pow( normalizing_radius, 3. );
        for (const auto& [key, value] : PARAMS::PLANETARY_GM){
            PARAMS::PLANETARY_GM[key] /= pow( normalizing_radius, 3.0 );
        }
    }

	void normalizeVariables(const double normalizing_radius) {
		// Normalize the distance
    	PARAMS::RP /= normalizing_radius;
		PARAMS::RA /= normalizing_radius;
    	PARAMS::SMAP /= normalizing_radius;
    	PARAMS::R /= normalizing_radius;
    	PARAMS::RS /= normalizing_radius;
    	// Gravitational parameter
    	PARAMS::MUP /= pow( normalizing_radius, 3. );
    	PARAMS::targetGM /= pow( normalizing_radius, 3. );
    	PARAMS::hostGM /= pow( normalizing_radius, 3. );
		for (const auto& [key, value] : PARAMS::PLANETARY_GM){
		    PARAMS::PLANETARY_GM[key] /= pow( normalizing_radius, 3.0 );
		}
    }

	void restoreVariables() {
		PARAMS::RA = _RA;
		PARAMS::RP = _RP;
		PARAMS::SMAP = _SMAP;
		PARAMS::RS = _RS;
		PARAMS::MUP = _MUP;
		PARAMS::targetGM = _targetGM;	
		PARAMS::hostGM = _hostGM;
		PARAMS::R = _R;
		PARAMS::PLANETARY_GM = _PLANETARY_GM;
	}

} // PARAMS:: namespace.

#endif
