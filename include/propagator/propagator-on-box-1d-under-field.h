#ifndef _PROPAGATOR_ON_BOX_1D_UNDER_FIELD_H_
#define _PROPAGATOR_ON_BOX_1D_UNDER_FIELD_H_

#include <complex>

#include "propagator-on-box-1d.h"


class Propagator_on_Box_1D_under_field : public Propagator_on_Box_1D {

private:
	double (*pAfunc)(double);
	double charge;
	double *M1 = NULL, *D1 = NULL;
	double *M1UA_half_dt_forward, *M1UA_half_dt_backward;

public:
	Propagator_on_Box_1D_under_field(
			size_t Nx, double dx, double *Vx, double (*pAfunc)(double), 
			double hbar=1., double mass=1., double charge=-1.);
	~Propagator_on_Box_1D_under_field();
	int propagate_under_field(std::complex<double> *wf, double dt, double t);
};

#endif // _PROPAGATOR_ON_BOX_1D_UNDER_FIELD_H_
