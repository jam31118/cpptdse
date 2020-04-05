import numpy as np
from numpy import exp

from qprop.param.file import ParameterFile
param = ParameterFile("./in.param")
Nx, dx = (param[key] for key in ("Nx", "dx"))

#Nx_tot = 1 + Nx + 1
xarr = 0. + dx * np.arange(1,Nx+1)
xmin, xmax = xarr[[0,-1]]
xmid = 0.2 * (xmin + xmax)
wf_t0 = exp(-(xarr-xmid)**2).astype(np.complex)
kx = -5.
wf_t0 *= exp(1.j*kx*xarr)
wf_t0[[0,-1]] = 0.

from tdse.propagator.box1d import Wavefunction_Uniform_1D_Box
Wavefunction_Uniform_1D_Box.normalize(wf_t0, dx)

wf_t0.tofile("wf-t0.bin")
