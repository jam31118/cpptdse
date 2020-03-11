#include <complex>
#include <cstdlib>
#include <iostream>

#include "param.h"

#include "../../../include/propagator/1d.hh"

int main() {
	
	ParamFile param;
	try { param = ParamFile("in.param"); }
	catch (std::exception& e) { 
		std::cerr << "[ERROR] Failed to construct parameter file object\n";
		std::cerr << "[ERROR] Error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch (...) { 
		std::cerr << "[ERROR] Failed during default constructor of ParamFile\n"; 
		return EXIT_FAILURE;
	}

	const size_t Nx = param.get_int("Nx");
	const double dx = param.get_double("dx");
	const double dt = param.get_double("dt");

	double *Vx = new double[Nx];
	for (double *pVx=Vx, *pVx_max=Vx+Nx; pVx<pVx_max; ++pVx) { *pVx = 0.; }

	Propagator_on_Box_1D prop = Propagator_on_Box_1D(Nx, dx, Vx);

	std::complex<double> wf[Nx];

	prop.propagate(wf, dt, 1);
	

	// Release resources
	prop.~Propagator_on_Box_1D();
	delete [] Vx;

	// Quit main program
	return EXIT_SUCCESS;
}

