#include <complex>
#include <cstdlib>
#include <iostream>
#include <random>
#include <fstream>
#include <algorithm>
#include <cmath>

#include "param.h"

#include "../../../include/propagator/1d.hh"
#include "../../../include/array.hh"

class Wavefunction_on_Box_1D {
	size_t Nx;
	double dx;
public:
	Wavefunction_on_Box_1D(size_t Nx, double dx);
	static double norm_sq(std::complex<double> *wf, size_t Nx, double dx);
	double norm_sq(std::complex<double> *wf);
	static int normalize(std::complex<double> *wf, size_t Nx, double dx);
	int normalize(std::complex<double> *wf);
};

Wavefunction_on_Box_1D::Wavefunction_on_Box_1D(size_t Nx, double dx): 
	Nx(Nx), dx(dx) {
	if (!(Nx > 0)) { throw "`Nx` should be positive"; }
	if (!(dx > 0)) { throw "`dx` should be real"; }
};

double Wavefunction_on_Box_1D::norm_sq(std::complex<double> *wf) {
	return this->norm_sq(wf, this->Nx, this->dx);
//	double _norm_sq = 0.0;
//	std::complex<double> _wf;
//	for (std::complex<double> *pwf=wf, *pwfmax=wf+Nx; pwf<pwfmax; ++pwf)
//	{ 
//		_wf = *pwf;
//		_norm_sq += std::real( std::conj(_wf) * _wf ); 
//	} _norm_sq *= dx;
//	return _norm_sq;
}

double Wavefunction_on_Box_1D::norm_sq(
		std::complex<double> *wf, size_t Nx, double dx) {
	double _norm_sq = 0.0;
	std::complex<double> _wf;
	for (std::complex<double> *pwf=wf, *pwfmax=wf+Nx; pwf<pwfmax; ++pwf)
	{ 
		_wf = *pwf;
		_norm_sq += std::real( std::conj(_wf) * _wf ); 
	} _norm_sq *= dx;
	return _norm_sq;
}

template <typename T1, typename T2>
int array_mul_scalar(T1 *a, size_t N, T2 c) {
	for (T1 *pa=a, *pamax=a+N; pa<pamax; ++pa) { *pa *= c; }
	return EXIT_SUCCESS;
}


int Wavefunction_on_Box_1D::normalize(
		std::complex<double> *wf, size_t Nx, double dx) {
	double _norm_sq = Wavefunction_on_Box_1D::norm_sq(wf, Nx, dx);
	return array_mul_scalar(wf, Nx, 1./std::sqrt(_norm_sq));	
}

int Wavefunction_on_Box_1D::normalize(std::complex<double> *wf) {
	return this->normalize(wf, this->Nx, this->dx);
}



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
	const size_t Nt = param.get_int("Nt");

	double *Vx = new double[Nx];
	for (double *pVx=Vx, *pVx_max=Vx+Nx; pVx<pVx_max; ++pVx) { *pVx = 0.; }

	Propagator_on_Box_1D prop = Propagator_on_Box_1D(Nx, dx, Vx);

	std::complex<double> *wf = new std::complex<double>[Nx];

	// Assign to random numbers
//	std::random_device rd;
//	std::mt19937 gen(rd());
//	std::uniform_real_distribution<double> unidis(-1., 1.);
//	for (std::complex<double> *pwf=wf, *pwfmax=wf+Nx; pwf<pwfmax; ++pwf)
//	{ *pwf = unidis(gen); } // imaginary part is set to zero by default


	const size_t Nx_tot = 1+Nx+1;
	const double L = dx * (Nx_tot-1);
	for (size_t i=0; i<Nx; ++i) {
		wf[i] = std::sin(M_PI / L * (i+1)*dx);
	}	


	Wavefunction_on_Box_1D::normalize(wf, Nx, dx);


	// Prepare storage for time-dependent wavefunction
	std::complex<double> *wf_t_1d = new std::complex<double>[Nx*Nt];
	std::complex<double> **wf_t = new std::complex<double>*[Nt];
	for (
			std::complex<double> **pwf_t=wf_t,**pwf_t_max=wf_t+Nt,*pwf_t_1d=wf_t_1d; 
			pwf_t<pwf_t_max; 
			++pwf_t, pwf_t_1d += Nx) 
	{ *pwf_t = pwf_t_1d; }
	
	std::complex<double> *wf_max = wf + Nx;
	std::copy(wf, wf_max, wf_t[0]);

	for (size_t it=0; it<Nt-1; ++it) {
		prop.propagate(wf, dt, 1);

//		std::cout << "wf[" << it+1 << "]= \n";
//		print_array(wf, Nx);

		std::copy(wf, wf_max, wf_t[it+1]);
	}
	

	// Write to output file
	std::string wf_t_file_name("wf_t.bin");
	std::ofstream wf_t_file(wf_t_file_name, std::ios::binary);
	if (!wf_t_file.is_open()) {
		std::cerr << "[ERROR] Failed to open file for `wf_t`\n";
		return EXIT_FAILURE;
	}
	wf_t_file.write( (char *) wf_t_1d, Nt*Nx*sizeof(std::complex<double>));
	wf_t_file.close();	
	std::cout << "[ LOG ] Wavefunction file written to: " 
		<< wf_t_file_name << std::endl; 


	// Release resources
	prop.~Propagator_on_Box_1D();
	delete [] Vx;
	delete [] wf;

	delete [] wf_t;
	delete [] wf_t_1d;

	// Quit main program
	return EXIT_SUCCESS;
}

