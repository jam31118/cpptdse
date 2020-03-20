#include <complex>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "param.h"

#include "../../../include/propagator/propagator-on-box-1d.h"
#include "../../../include/wf/wavefunction-on-box-1d.h"
#include "../../../include/array.h"


#ifdef FIELD

#include "../../../include/propagator/propagator-on-box-1d-under-field.h"

//double A_func(double) { return 0.; }

double A_func(double t) {
	double A0 = 0.029638422288250095, w = 0.05695419063519442, nc = 2;
	double ww = w / (2.*nc);
	double duration = M_PI/ww;
	double start_time = 0.0;
	double end_time = start_time + duration;
	if (t > end_time) { return 0.; }
	else { return A0 * sin(ww*t) * sin(w*t); }
}

#endif // FIELD


int main() {
	
	// Extract parameters from file
	//
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


	// Construct a time-independent potential array
	double *Vx = new double[Nx];
	set_to_zeros(Vx, Nx);

	// Construct a propagator object
#ifdef FIELD
	Propagator_on_Box_1D_under_field prop(Nx, dx, Vx);
#else // FIELD
	Propagator_on_Box_1D prop(Nx, dx, Vx);
#endif // FIELD

	
	// Prepare initial state
	//
	std::complex<double> *wf = new std::complex<double>[Nx];
	set_to_randoms(wf, Nx); // Initialize wavefunction array to random numbers
	// Propagate to lowest energy possible from given wf
	if (prop.propagate_to_ground_state(wf, dt, 5000, 1e-13) != EXIT_SUCCESS) {
		std::cerr << "[ERROR] Failed to propagate to ground state\n";
		return EXIT_FAILURE;		
	} std::cout << "[ LOG ] COMPLETE: A propagation to ground state\n";


	// Prepare storage for time-dependent wavefunction
	// 
	std::complex<double> *wf_t_1d = new std::complex<double>[Nx*Nt];
	std::complex<double> **wf_t = new std::complex<double>*[Nt];
	set_2d_view_of_1d(wf_t, wf_t_1d, Nt, Nx);
	

	// Main Propagation
	// 
#ifdef FIELD
	double t = 0.;
	double *At = new double[Nt];
	At[0] = A_func(t);
#endif // FIELD
	std::complex<double> *wf_max = wf + Nx;
	std::copy(wf, wf_max, wf_t[0]);
	for (size_t it=0; it<Nt-1; ++it) {
#ifdef FIELD
		prop.propagate_under_field(wf, dt, A_func(t+0.5*dt)); // should be t+0.5*dt
		t += dt;
		At[it+1] = A_func(t);
#else // FIELD
		prop.propagate(wf, dt, 1);
#endif // FIELD
//		std::cout << "wf[" << it+1 << "]= \n";
//		print_array(wf, Nx);
		std::copy(wf, wf_max, wf_t[it+1]);
	} std::cout << "[ LOG ] COMPLETE: A propagation with real timestep\n";
	

	// Write to output file
	//
	std::string wf_t_file_name("wf_t.bin");
	std::ofstream wf_t_file(wf_t_file_name, std::ios::binary);
	if (!wf_t_file.is_open()) {
		std::cerr << "[ERROR] Failed to open file for `wf_t`\n";
		return EXIT_FAILURE;
	}
	wf_t_file.write( (char *) wf_t_1d, Nt*Nx*sizeof(std::complex<double>));
	wf_t_file.close();	
	std::cout << "[ LOG ] COMPLETE: Writing wavefunction to a file: " 
		<< wf_t_file_name << std::endl; 

#ifdef FIELD
	std::ofstream vecpot_t_file("vecpot-t.bin");
	vecpot_t_file.write((char *) At, Nt*sizeof(double));
	vecpot_t_file.close();
#endif // FIELD




	// Release resources
	//
#ifdef FIELD
	delete [] At;	
#endif // FIELD
	delete [] Vx;
	delete [] wf;

	delete [] wf_t;
	delete [] wf_t_1d;


	// Quit main program
	return EXIT_SUCCESS;
}

