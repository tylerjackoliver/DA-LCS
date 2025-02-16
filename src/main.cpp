#include <iostream>                 // I/O
#include <mpi.h>
#include <dace/dace.h>              // Differential algebra
#include <boost/numeric/odeint.hpp> // Numerical integration
#include "Params.hpp"               // Simulation parameters
#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <fstream>
#include "eigenvectors.h"
#include <chrono>
#include <random>
#include "helicityIntegrator.h"
#include "drivers.h"

int main(int argv, char* argc[])
{
    int poolSize, myRank;
    MPI_Init(nullptr, nullptr);
    DACE::DA::init(2, 3);
    DACE::DA::setEps(0.0);
    PARAMS::initialiseVariables();
	PARAMS::normalizeVariables(PARAMS::RP);

	/* Driver code goes here. */

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
