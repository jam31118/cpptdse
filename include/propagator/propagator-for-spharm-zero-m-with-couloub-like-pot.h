#ifndef _PROPAGATOR_FOR_SPHARM_ZERO_M_WITH_COULOMB_LIKE_POT_H_
#define _PROPAGATOR_FOR_SPHARM_ZERO_M_WITH_COULOMB_LIKE_POT_H_

#include "propagator-for-spharm-zero-m.h"

class Propagator_for_spharm_zero_m_with_coulomb_like_pot 
: public Propagator_for_spharm_zero_m 
{
public:
	Propagator_for_spharm_zero_m_with_coulomb_like_pot(
			size_t Nr, double dr, size_t Nl, double *Vr, double Z, 
			double hbar=1., double mass=1.);
};

#endif // _PROPAGATOR_FOR_SPHARM_ZERO_M_WITH_COULOMB_LIKE_POT_H_
