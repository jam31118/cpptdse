#include "../../include/propagator/\
propagator-for-spharm-zero-m-with-couloub-like-pot.h"


Propagator_for_spharm_zero_m_with_coulomb_like_pot
::Propagator_for_spharm_zero_m_with_coulomb_like_pot(
		size_t Nr, double dr, size_t Nl, double *Vr, double Z, 
		double hbar, double mass): 
	Propagator_for_spharm_zero_m(Nr, dr, Nl, Vr, hbar, mass)
{
	//// NOTE on correction for coulomb-like potential in a form of V(r)=-Z/r	
	//
	// For l=0 :
	// - ReVl = ReVr
	// - M2ReH[il==0] = Kr(corrected) + M2*ReVl
	//                = Kr(corrected) + M2*ReVr
	//                = (-hbar**2/(2*mass)) * Dr(corrected) + M2*ReVr
	// - M2ReH[il==0][diag][0] = (-hbar**2/(2*mass)) * Dr(corrected)[diag][0] 
	//                           + M2[diag][0] * ReVr[0]
	// - Dr(corrected)[diag][0] = (delta)
	// - (delta) \equiv -2/dr**2 * (1 - Z * dr / (12 - 10 * Z * dr))
	const double delta = -2./(dr*dr) * (1 - Z * dr / (12. - 10. * Z * dr));
	const size_t il = 0;  // index for zero l
	const size_t ir = 0;
	// index of diagonal line (length of Nr) in tridiag array (length of 3*Nr), 
	// which is arranged in order of lower-offdiag, diag, upper-offdiag.
	const size_t i_diag = Nr + ir;
	M2ReH[il][i_diag] = (-hbar*hbar/(2*mass)) * delta + M2[i_diag] * Vr[ir];
}
