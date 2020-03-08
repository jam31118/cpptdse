#include <complex>
#include <cstdlib>

#include "../../../include/propagator/1d.hh"

int main() {
	const size_t Nx = 9;
	double Vx[Nx] = {0.0};
	Propagator_on_Box_1D prop = Propagator_on_Box_1D(9, 0.2, Vx);	

	std::complex<double> wf[Nx];

	prop.propagate(wf, 0.05, 1);
	
	prop.~Propagator_on_Box_1D();

	return EXIT_SUCCESS;
}
