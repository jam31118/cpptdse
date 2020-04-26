#ifndef _PROPAGATOR_FOR_SPHARM_ZERO_M_H_
#define _PROPAGATOR_FOR_SPHARM_ZERO_M_H_

#include <complex>
#include "../wf/wavefunction-for-spharm-zero-m.h"

class Propagator_for_spharm_zero_m {

	size_t Nr;
	double dr;
	size_t Nl;
	const double *Vr;
	double hbar, mass;

	size_t N_tridiag;
	double *M2ReH_1d;
	std::complex<double> *U_1d, **U, *Uinv_1d, **Uinv;

protected:
	double *M2;
	double **M2ReH;

public:

	Propagator_for_spharm_zero_m(
		const size_t Nr, const double dr, const size_t Nl, const double *Vr, 
		const double hbar=1., const double mass=1.);
	~Propagator_for_spharm_zero_m();
	int eval_time_evol_unitary_for_real_timestep(double dt);
	int eval_time_evol_unitary_for_imag_timestep(
			double imag_dt, size_t nonzero_l_num);
	int propagate(
			std::complex<double> *const wf, const double dt, const size_t Nt=1);
	int propagate_to_ground_state(
		std::complex<double> *const wf, const double dt,
		const size_t Nt_max, const double stop_thres, const size_t nonzero_l_num);

	Wavefunction_for_spharm_zero_m *wf;
};


#endif // _PROPAGATOR_FOR_SPHARM_ZERO_M_H_
