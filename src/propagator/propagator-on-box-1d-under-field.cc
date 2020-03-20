#include "../../include/propagator/propagator-on-box-1d-under-field.h"

#include "../../include/tridiag.h"
#include "../../include/array.h"


Propagator_on_Box_1D_under_field::Propagator_on_Box_1D_under_field(
		size_t Nx, double dx, double *Vx, double (*pAfunc)(double),
		double hbar, double mass, double charge) : 
	Propagator_on_Box_1D(Nx, dx, Vx, hbar, mass), pAfunc(pAfunc), charge(charge)
{

	M1 = new double[N_tridiag];
	eval_M1_tridiag(M1, Nx);

	D1 = new double[N_tridiag];
	if (EXIT_SUCCESS != eval_D1_tridiag(D1, dx, Nx))
	{ throw "Failed to evaluate D1 tridiagonal"; }

	M1UA_half_dt_forward = new double[N_tridiag]; 
	M1UA_half_dt_backward = new double[N_tridiag]; 

}


Propagator_on_Box_1D_under_field::~Propagator_on_Box_1D_under_field() {
	delete [] M1;
	delete [] D1;
	delete [] M1UA_half_dt_forward;
	delete [] M1UA_half_dt_backward;
}


int Propagator_on_Box_1D_under_field::propagate_under_field(
		std::complex<double> *wf, double dt, double t) {

	//// Prepare time evolution unitary tridiagonals
	//
	// For H0: the Hamiltonian without field
	// The following method evaluates `U_forward` and `U_backward`
	eval_time_evol_unitary_for_real_timestep(0.5*dt);
	std::complex<double> 
		*M2U0_quarter_dt_forward = U_forward, // aliasing 
		*M2U0_quarter_dt_backward = U_backward; // aliasing
	//
	// For HA: the Hamiltonian under electromagnatic field in velocity gauge
	const double At_at_half = (*pAfunc)(t+0.5*dt);
	const double coef = 0.5 * dt * charge / mass * At_at_half;
	v1_add_c_mul_v2(M1UA_half_dt_forward, M1, coef, D1, N_tridiag);
	v1_add_c_mul_v2(M1UA_half_dt_backward, M1, -coef, D1, N_tridiag);

	//// Propagation using split-operator method
	//
	std::complex<double> *wf_mid = new std::complex<double>[Nx];

	tridiag_forward(M2U0_quarter_dt_forward, wf, wf_mid, Nx);
	tridiag_backward(M2U0_quarter_dt_backward, wf, wf_mid, Nx);

	tridiag_forward(M1UA_half_dt_forward, wf, wf_mid, Nx);
	tridiag_backward(M1UA_half_dt_backward, wf, wf_mid, Nx);

	tridiag_forward(M2U0_quarter_dt_forward, wf, wf_mid, Nx);
	tridiag_backward(M2U0_quarter_dt_backward, wf, wf_mid, Nx);

	delete [] wf_mid;

	return EXIT_SUCCESS;
}

