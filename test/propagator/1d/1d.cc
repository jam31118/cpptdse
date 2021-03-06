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




#ifdef IMAGPOT

#include <cmath>
int eval_imagpot_powered(
		double *imagpot, size_t Nx, size_t Nwidth, double V0) {
	if (V0 < 0) { return EXIT_FAILURE; }
	if (Nx < (2*Nwidth+1)) { return EXIT_FAILURE; }
	double *pimagpot = imagpot;
	size_t n = 0; size_t ord = 7;

	double double_Nwidth = Nwidth;
	double double_grid_num;
//		std::cout << ((double) Nwidth) << std::endl;
//		std::cout <<  double_Nwidth << std::endl;
	for ( ; n < Nwidth; ++n, ++pimagpot)
	{ 
		double_grid_num = (int)(n - (Nwidth - 1));
//		std::cout << double_grid_num << double_Nwidth << std::endl;
		*pimagpot = V0 * pow(double_grid_num / double_Nwidth, ord); 
	}
	for ( ; n < Nx - Nwidth; ++n, ++pimagpot)
	{ *pimagpot = 0.; }
	for ( ; n < Nx; ++n, ++pimagpot)
	{ 
		double_grid_num = (int)(Nx - Nwidth - n);
		*pimagpot = V0 * pow(double_grid_num / double_Nwidth, ord); 
	}
	return EXIT_SUCCESS;
}

#endif // IMAGPOT




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


#ifdef IMAGPOT
	double *Vx_imag = new double[Nx];
//	set_to_zeros(Vx_imag, Nx);
	const size_t N_imagpot_width = param.get_int("imagpot-grid-num");
	const double imagpot_ampl = param.get_double("imagpot-ampl");
	if (EXIT_SUCCESS != eval_imagpot_powered(
				Vx_imag, Nx, N_imagpot_width, imagpot_ampl)) {
		std::cerr << "[ERROR] Failed to evaluate imaginary scalar potential\n";
		return EXIT_FAILURE;
	}
#endif // IMAGPOT


	// Construct a propagator object
#ifdef FIELD
	Propagator_on_Box_1D_under_field prop(Nx, dx, Vx, &A_func);
#elif IMAGPOT
	Propagator_on_Box_1D_with_imag_pot prop(Nx, dx, Vx, Vx_imag);
#else
	Propagator_on_Box_1D prop(Nx, dx, Vx);
#endif // FIELD

	
	// Prepare initial state
	//
	std::complex<double> *wf = new std::complex<double>[Nx];

	// Try to get the name of file for initial wavefunction data, if any.
	std::string wf_t0_fname;
	try { wf_t0_fname = param.get_string("wf-t0-file"); }
	catch (...) { 
		std::cout << "[ LOG ] No file for initial wavefunction found. "
			"Falling back to search for ground state\n";
	}
	if (!wf_t0_fname.empty()) {
		try { 
			std::ifstream wf_t0_file(wf_t0_fname, std::ios::binary);
			wf_t0_file.read((char *) wf, Nx * sizeof(std::complex<double>));
			wf_t0_file.close();
		} catch (...) {
			std::cerr << "[ERROR] Failed to read initial wavefunction from: " 
				<< wf_t0_fname << std::endl;
			return EXIT_FAILURE;
		}
	} else {
		set_to_randoms(wf, Nx); // Initialize wavefunction array to random numbers
		// Propagate to lowest energy possible from given wf
		if (prop.propagate_to_ground_state(wf, dt, 20000, 1e-13) != EXIT_SUCCESS) {
			std::cerr << "[ERROR] Failed to propagate to ground state\n";
			return EXIT_FAILURE;		
		} std::cout << "[ LOG ] COMPLETE: A propagation to ground state\n";
	}


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
		prop.propagate_under_field(wf, dt, t);
		t += dt;
		At[it+1] = A_func(t);
#else // FIELD
		prop.propagate(wf, dt, 1);
#endif // FIELD
//		std::cout << "wf[" << it+1 << "]= \n";
//		print_array(wf, Nx);
		std::copy(wf, wf_max, wf_t[it+1]);
	} std::cout << "[ LOG ] COMPLETE: A propagation with real timestep\n";
	



	//// Write to output file
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


	std::ofstream Vx_file("scalarpot_real.bin");
	Vx_file.write((char *) Vx, Nx*sizeof(double));
	Vx_file.close();

#ifdef IMAGPOT
	std::ofstream Vx_imag_file("scalarpot_imag.bin");
	Vx_imag_file.write((char *) Vx_imag, Nx*sizeof(double));
	Vx_imag_file.close();
#endif // IMAGPOT

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
#ifdef IMAGPOT
	delete [] Vx_imag;
#endif // IMAGPOT
	delete [] wf;

	delete [] wf_t;
	delete [] wf_t_1d;


	// Quit main program
	return EXIT_SUCCESS;
}

