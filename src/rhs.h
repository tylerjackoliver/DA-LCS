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

#ifndef ACROBAT_RHS_H
#define ACROBAT_RHS_H

#include <dace/dace_s.h>
#include "Params.hpp"
#include <vector>
#include <fstream>
extern "C"
{
#include <SpiceUsr.h>
}

namespace RHS
{

    void forceFunction(const DACE::AlgebraicVector<DACE::DA>& x0, DACE::AlgebraicVector<DACE::DA>& x, double t)
    {
        if (x.size() != x0.size()) x.resize( x0.size() );
        x[0] = x0[3];
        x[1] = x0[4];
        x[2] = x0[5];
        double fR = 1. / ( 1 + PARAMS::EP * cos( PARAMS::T0 + t ) );
        DACE::DA r1 = sqrt( (x0[0] + PARAMS::MU) * (x0[0] + PARAMS::MU) + x0[1]*x0[1] + x0[2] * x0[2] );
        DACE::DA r2 = sqrt( (x0[0] + PARAMS::MU - 1) * (x0[0] + PARAMS::MU - 1) + x0[1]*x0[1] + x0[2]*x0[2] );
        DACE::DA r13 = r1 * r1 * r1;
        DACE::DA r23 = r2 * r2 * r2;
        DACE::DA aX = 2.0 * x0[4] + fR * ( x0[0] - (1 - PARAMS::MU) * (x0[0] + PARAMS::MU) * (1/r13) -
                                           PARAMS::MU * (x0[0] + PARAMS::MU - 1) * (1/r23) );
        DACE::DA aY = -2.0 * x0[3] + fR * ( x0[1] - (1 - PARAMS::MU) * x0[1] * (1/r13) - PARAMS::MU * x0[1] * (1/r23) );
        DACE::DA aZ = -x0[2] + fR * (x0[2] - (1 - PARAMS::MU) * x0[2] * (1/r13) - (PARAMS::MU * x0[2]) * (1/r23) );
        x[3] = aX;
        x[4] = aY;
        x[5] = aZ;
    }

    void forceFunctionObserver(const DACE::AlgebraicVector<DACE::DA>& x, double t)
    {
        ;
    }

    void forceFunction(const DACE::AlgebraicVector<double>& x0, DACE::AlgebraicVector<double>& x, double t)
    {
        if (x.size() != x0.size()) x.resize( x0.size() );
        x[0] = x0[3];
        x[1] = x0[4];
        x[2] = x0[5];
        double fR = 1. / ( 1 + PARAMS::EP * cos( PARAMS::T0 + t ) );
        double r1 = sqrt( (x0[0] + PARAMS::MU) * (x0[0] + PARAMS::MU) + x0[1]*x0[1] + x0[2] * x0[2] );
        double r2 = sqrt( (x0[0] + PARAMS::MU - 1) * (x0[0] + PARAMS::MU - 1) + x0[1]*x0[1] + x0[2]*x0[2] );
        double r13 = r1 * r1 * r1;
        double r23 = r2 * r2 * r2;
        double aX = 2.0 * x0[4] + fR * ( x0[0] - (1 - PARAMS::MU) * (x0[0] + PARAMS::MU) * (1/r13) -
                                           PARAMS::MU * (x0[0] + PARAMS::MU - 1) * (1/r23) );
        double aY = -2.0 * x0[3] + fR * ( x0[1] - (1 - PARAMS::MU) * x0[1] * (1/r13) - PARAMS::MU * x0[1] * (1/r23) );
        double aZ = -x0[2] + fR * (x0[2] - (1 - PARAMS::MU) * x0[2] * (1/r13) - (PARAMS::MU * x0[2]) * (1/r23) );
        x[3] = aX;
        x[4] = aY;
        x[5] = aZ;
    }

    void forceFunctionObserver(const DACE::AlgebraicVector<double>& x, double t)
    {
        ;
    }

}

#endif //ACROBAT_RHS_H
