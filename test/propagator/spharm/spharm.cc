#include <iostream>
#include <complex>
#include <fstream>

#include "param.h"

#include "../../../include/array.h"
#include "../../../include/propagator/propagator-for-spharm-zero-m.h"



#include <cmath>
int eval_ground_state_of_sphbox_in_spharm_basis(
		std::complex<double> *wf, size_t Nr, double dr, size_t Nl) 
{
	if (Nr < 2 || Nl < 1 || dr <= 0.) { return EXIT_FAILURE; }
	const double rmax = (Nr+1) * dr;
	const double norm_const = std::sqrt(2./rmax);
	const double dphi = M_PI*dr/rmax;
	double phi=0;
	std::complex<double> *pwfmax = wf + Nr;
	for (std::complex<double> *pwf=wf; pwf<pwfmax; ++pwf) {
		phi += dphi;
		*pwf = norm_const * std::sin(phi);
	}
	set_to_zeros(wf+Nr*1, Nr*(Nl-1));
	return EXIT_SUCCESS;
}



int main() {

	// Extract parameters from file
	//
	ParamFile *param = NULL;
	try { param = new ParamFile("in.param"); }
	catch (...) { 
		std::cerr << "[ERROR] Failed to construct param objects\n"; 
		return EXIT_FAILURE;
	}
	const size_t Nr = param->get_long("Nr");
	const double dr = param->get_double("dr");
	const size_t Nl = param->get_long("Nl");
	const size_t Nt = param->get_long("Nt");
	const double dt = param->get_double("dt");



	double *Vr = new double[Nr];
	set_to_zeros(Vr, Nr);

	Propagator_for_spharm_zero_m *prop = NULL;
	try { prop = new Propagator_for_spharm_zero_m(Nr, dr, Nl, Vr); }
	catch (const char *mesg) {
		std::cerr << "[ERROR] Failed to construct propagator with message: \n";
		std::cerr << "[ERROR] '" << mesg << "'\n";
		return EXIT_FAILURE;
	}

	
	// Prepare initial state
	//
	const size_t wf_length = prop->wf->length();
	std::complex<double> *wf = new std::complex<double>[wf_length];
	std::complex<double> *wf_max = wf + wf_length;
	
	int stat = eval_ground_state_of_sphbox_in_spharm_basis(wf, Nr, dr, Nl);
	if (EXIT_SUCCESS != stat) {
		std::cerr << "[ERROR] Failed to evaluate ground state\n";
		return EXIT_FAILURE;
	}


	// Copy the initial wavefunction into a separate array
	std::complex<double> *wf_t0 = new std::complex<double>[wf_length];
	std::copy(wf, wf_max, wf_t0);



	// Prepare storage for time-dependent wavefunction
	// 
	std::complex<double> *wf_t_1d = new std::complex<double>[Nt*wf_length];
	std::complex<double> **wf_t = new std::complex<double>*[Nt];
	set_2d_view_of_1d(wf_t, wf_t_1d, Nt, Nr);


	//// Propagate wavefunction
	//
	double t = 0.;
	double *tarr = new double[Nt];
	tarr[0] = t;
	std::copy(wf, wf_max, wf_t[0]);
	for (size_t it=0; it<Nt-1; ++it) {
		prop->propagate(wf, dt, 1);
		t += dt;
		tarr[it+1] = t;
		std::copy(wf, wf_max, wf_t[it+1]);
	} std::cout << "[ LOG ] COMPLETE: A propagation with real timestep\n";


	//// Save data into files
	//
	// Save radial coordinate array
	double *rarr = new double[Nr];
	prop->wf->eval_radial_arr(rarr);
	std::ofstream rarr_file("rarr.bin", std::ios::binary);
	rarr_file.write((char *) rarr, sizeof(double) * Nr);
	rarr_file.close();
	delete [] rarr;

	// Save the initial wavefunction
	std::ofstream wf_t0_file("wf-t0.bin", std::ios::binary);
	wf_t0_file.write((char *) wf_t0, sizeof(std::complex<double>) * wf_length);
	wf_t0_file.close();

	// Save the static potential
	std::ofstream Vr_file("Vr_real.bin");
	Vr_file.write((char *) Vr, Nr*sizeof(double));
	Vr_file.close();

	// Save the time-dependent wavefunction
	std::string wf_t_file_name("wf_t.bin");
	std::ofstream wf_t_file(wf_t_file_name, std::ios::binary);
	wf_t_file.write((char *) wf_t_1d, Nt*wf_length*sizeof(std::complex<double>));
	wf_t_file.close();	
	std::cout << "[ LOG ] COMPLETE: Writing wavefunction to a file: " 
		<< wf_t_file_name << std::endl; 

	// Save the time array
	std::ofstream tarr_file("tarr.bin", std::ios::binary);
	tarr_file.write((char *) tarr, sizeof(double) * Nt);
	tarr_file.close();
	delete [] tarr;


	delete [] wf_t_1d;
	delete [] wf_t;
	delete [] wf_t0;
	delete [] wf;
	delete [] Vr;
	delete prop;
	delete param;

	return EXIT_SUCCESS;
}
