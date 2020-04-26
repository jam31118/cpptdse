#include "../../include/propagator/propagator-for-spharm-zero-m.h"

#include "../../include/tridiag.h"
#include "../../include/array.h"
#include "../../include/wf/wavefunction-for-spharm-zero-m.h"


Propagator_for_spharm_zero_m::Propagator_for_spharm_zero_m(
		size_t Nr, double dr, size_t Nl, double *Vr, double hbar, double mass): 
	Nr(Nr), dr(dr), Nl(Nl), Vr(Vr), hbar(hbar), mass(mass)
{

	if (dr <= 0.) { throw "Negative `dr` given"; }
	if (Nr < 2) { throw "`Nr` should be two or bigger"; }
	if (Nl < 1) { throw "`Nl` should be one or bigger"; }

	wf = new Wavefunction_for_spharm_zero_m(Nr, dr, Nl);

	N_tridiag = 3 * Nr;

	// Allocate memory for time evolution unitary operators and set 2D view
	U_1d = new std::complex<double>[Nl*N_tridiag];
	U = new std::complex<double>*[Nl];
	set_2d_view_of_1d(U, U_1d, Nl, N_tridiag);
	Uinv_1d = new std::complex<double>[Nl*N_tridiag];
	Uinv = new std::complex<double>*[Nl];
	set_2d_view_of_1d(Uinv, Uinv_1d, Nl, N_tridiag);
	
	// Evaluate M2 matric for 'Numerov boost'
	M2 = new double[N_tridiag];
	if (EXIT_SUCCESS != eval_M2_tridiag(M2, Nr)) 
	{ throw "Failed to evaluate M2 tridiag"; }


	//// Construct M2ReH matrix for constructing time evolution matrix
	//
	// Allocate memory for M2ReH matrix and set 2D view
	M2ReH_1d = new double[Nl*N_tridiag];
	M2ReH = new double*[Nl];
	set_2d_view_of_1d(M2ReH, M2ReH_1d, Nl, N_tridiag);
	// Constant
	const double hbar_sq_over_2mass = hbar*hbar/(2.*mass);
	// Construct Kr matrix
	double *const Kr = new double[N_tridiag]; // == -hbar**2/(2*mass)*D2
	eval_D2_tridiag(Kr, dr, Nr, -hbar_sq_over_2mass);
	// Storage for matrix Vl and M2ReVl
	double *Vl = new double[Nr];
	double *const M2ReVl = new double[N_tridiag];
	// Construct array: 1/r^2
	double *one_over_r_sq = new double[Nr];
	for (double *pr=one_over_r_sq, *prmax=one_over_r_sq+Nr, rval=dr;
			pr<prmax; ++pr, rval+=dr) { *pr = 1./(rval*rval); }
	// Evaluate ReVl and M2ReVl for M2ReH = Kr + M2*ReVl
	// where ReVl = ReVr + (hbar / (2mass)) * l * (l+1) / r^2
	for (size_t il=0; il<Nl; ++il) {
		const size_t l = il;  // should be modified for nonzero m
		v1_add_c_mul_v2(Vr, one_over_r_sq, hbar_sq_over_2mass*l*(l+1), Vl, Nr);	
		tridiag_mul_diag(M2, Vl, M2ReVl, Nr);
		add(Kr, M2ReVl, M2ReH[il], N_tridiag);
	}
	delete [] one_over_r_sq;
	delete [] M2ReVl;
	delete [] Vl;
	delete [] Kr;
}




Propagator_for_spharm_zero_m::~Propagator_for_spharm_zero_m() {

	delete [] U;
	delete [] U_1d;
	delete [] Uinv;
	delete [] Uinv_1d;

	delete [] M2;

	delete [] M2ReH;
	delete [] M2ReH_1d;

	delete wf;
}




int Propagator_for_spharm_zero_m
::eval_time_evol_unitary_for_real_timestep(double dt) {
		
	// Construct unitary propagators
	const double _c = - 0.5 * dt / hbar;
	double *const pM2max = M2 + N_tridiag;
	for (size_t il=0; il<Nl; ++il) 
	{
		double *pM2, *pM2ReHl;
		std::complex<double> *Ul=U[il], *Ul_inv=Uinv[il];
		std::complex<double> *pUl, *pUl_inv;
		for (pUl=Ul, pUl_inv=Ul_inv, pM2=M2, pM2ReHl=M2ReH[il];
				pM2 < pM2max;
				++pUl, ++pUl_inv, ++pM2, ++pM2ReHl) 
		{
			const double _M2 = *pM2;
			const double _c_M2ReHl = _c * (*pM2ReHl);
			*pUl = {_M2, _c_M2ReHl};
			*pUl_inv = {_M2, -_c_M2ReHl};
		}
		Ul[0] = 0., Ul_inv[0] = 0.;
		Ul[N_tridiag-1] = 0., Ul_inv[N_tridiag-1] = 0.;
	}
	return EXIT_SUCCESS;
}




int Propagator_for_spharm_zero_m
::eval_time_evol_unitary_for_imag_timestep(
		double imag_dt, size_t nonzero_l_num) 
{
		
	// Construct unitary propagators
	const double _c = - 0.5 * imag_dt / hbar;
	double *const pM2max = M2 + N_tridiag;
	for (size_t il=0; il<nonzero_l_num; ++il) 
	{
		double *pM2, *pM2ReHl;
		std::complex<double> *Ul=U[il], *Ul_inv=Uinv[il];
		std::complex<double> *pUl, *pUl_inv;
		for (pUl=Ul, pUl_inv=Ul_inv, pM2=M2, pM2ReHl=M2ReH[il];
				pM2 < pM2max;
				++pUl, ++pUl_inv, ++pM2, ++pM2ReHl) 
		{
			const double _M2 = *pM2;
			const double _c_M2ReHl = _c * (*pM2ReHl);
			*pUl = {_M2 + _c_M2ReHl, 0.};
			*pUl_inv = {_M2 -_c_M2ReHl, 0.};
		}
		Ul[0] = 0., Ul_inv[0] = 0.;
		Ul[N_tridiag-1] = 0., Ul_inv[N_tridiag-1] = 0.;
	}
	return EXIT_SUCCESS;
}



int Propagator_for_spharm_zero_m::propagate(
		std::complex<double> *const wf, const double dt, const size_t Nt) 
{
	// Construct unitary propagators
	eval_time_evol_unitary_for_real_timestep(dt);

	// Propagate
	std::complex<double> *wf_mid = new std::complex<double>[Nr];
	for (size_t it=0; it<Nt; ++it) {
		for (size_t il=0; il<Nl; ++il) {
			tridiag_forward(U[il], wf+il*Nl, wf_mid, Nr);	
			tridiag_backward(Uinv[il], wf+il*Nl, wf_mid, Nr);
		}
	}
	delete [] wf_mid;
	
	return EXIT_SUCCESS;
}



int Propagator_for_spharm_zero_m::propagate_to_ground_state(
		std::complex<double> *const wf, const double dt,
		const size_t Nt_max, const double stop_thres, const size_t nonzero_l_num) {
	
	const size_t nonzero_wf_length = nonzero_l_num * Nr;
	std::complex<double> *const nonzero_wf_max = wf + nonzero_wf_length;
	const size_t zero_l_num = Nl - nonzero_l_num;
	if (zero_l_num < 0) { return EXIT_FAILURE; }
	const size_t zero_wf_length = zero_l_num * Nr;

	// Set the rest of the wavefunction to zero
	set_to_zeros(wf + nonzero_wf_length, zero_wf_length);


	std::complex<double> *const nonzero_wf_prev 
		= new std::complex<double>[nonzero_wf_length];
	Wavefunction_for_spharm_zero_m::normalize(wf, Nr, dr, nonzero_l_num);
	std::copy(wf, nonzero_wf_max, nonzero_wf_prev);

	std::complex<double> *const nonzero_wf_diff
		= new std::complex<double>[nonzero_wf_length];

	double norm_diff;


	//// Propagate
	//
	eval_time_evol_unitary_for_imag_timestep(dt, nonzero_l_num);

	std::complex<double> *wf_mid = new std::complex<double>[Nr];
	size_t it;
	for (it=0; it<Nt_max; ++it) {

		// Propagate one step with an imaginary time step
		for (size_t il=0; il<nonzero_l_num; ++il) {
			tridiag_forward(U[il], wf+il*Nl, wf_mid, Nr);	
			tridiag_backward(Uinv[il], wf+il*Nl, wf_mid, Nr);
		}

		Wavefunction_for_spharm_zero_m::normalize(wf, Nr, dr, nonzero_l_num);
		substract(wf, nonzero_wf_prev, nonzero_wf_diff, nonzero_wf_length);
		norm_diff = Wavefunction_for_spharm_zero_m
			::norm_sq(nonzero_wf_diff, Nr, dr, nonzero_l_num);
		if (norm_diff < stop_thres) { break; }

		std::copy(wf, nonzero_wf_max, nonzero_wf_prev);
	}
	delete [] wf_mid;

	int return_code = EXIT_SUCCESS;
	if (it >= Nt_max) { 
		std::cerr << "[ERROR] Maximum iteration (=" << Nt_max << ") exceeded\n"; 
		// should not return from here 
		// since we need to proceed to complete resource release process
		return_code = EXIT_FAILURE;  
	}

	delete [] nonzero_wf_prev;
	delete [] nonzero_wf_diff;

	return EXIT_SUCCESS;
}


