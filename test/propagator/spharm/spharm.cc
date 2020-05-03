#include <iostream>
#include <complex>
#include <fstream>

#include "param.h"

#include "../../../include/array.h"

#ifdef COULOMB
#include "../../../include/propagator/\
propagator-for-spharm-zero-m-with-couloub-like-pot.h"
#else // COULOMB
#include "../../../include/propagator/propagator-for-spharm-zero-m.h"
#endif // COULOMB



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


	// Construct static potential
	double *Vr = new double[Nr];
#ifdef COULOMB
	const double Z = 1.;
	for (double r=dr, *pVr=Vr, *const pVr_max=Vr+Nr; pVr < pVr_max; ++pVr, r+=dr)
	{ *pVr = - Z / r; }
#else
	set_to_zeros(Vr, Nr);
#endif

	
	//// Construct a propagator
	//
#ifdef COULOMB
	Propagator_for_spharm_zero_m_with_coulomb_like_pot *prop = NULL;
	try { 
		prop = new Propagator_for_spharm_zero_m_with_coulomb_like_pot(
				Nr, dr, Nl, Vr, Z);	
	} catch (...) {
		std::cerr << "[ERROR] Failed to construct propagator\n";
		return EXIT_FAILURE;
	}
#else // COULOMB
	Propagator_for_spharm_zero_m *prop = NULL;
	try { prop = new Propagator_for_spharm_zero_m(Nr, dr, Nl, Vr); }
	catch (const char *mesg) {
		std::cerr << "[ERROR] Failed to construct propagator with message: \n";
		std::cerr << "[ERROR] '" << mesg << "'\n";
		return EXIT_FAILURE;
	}
#endif // COULOMB

	

	// Prepare initial state
	//
	const size_t wf_length = prop->wf->length();
	std::complex<double> *wf = new std::complex<double>[wf_length];
	std::complex<double> *wf_max = wf + wf_length;
	
	// Evaluate initial state by iteration : propagation with imaginary time
	set_to_randoms(wf, wf_length);
	prop->propagate_to_ground_state(wf, dt, 10000, 1e-14, 1);

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


	// Free memory
	//
	delete [] wf_t_1d;
	delete [] wf_t;
	delete [] wf_t0;
	delete [] wf;
	delete [] Vr;
	delete prop;
	delete param;

	return EXIT_SUCCESS;
}
