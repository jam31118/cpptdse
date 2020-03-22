import numpy as np
import matplotlib.pyplot as plt

from qprop.param.file import ParameterFile
param = ParameterFile("./in.param")
Nx = param['Nx']
Nt = param['Nt']

wf_t_1d = np.fromfile("./wf_t.bin", dtype=np.complex)
wf_t = wf_t_1d.reshape(Nt, Nx)
print(wf_t.shape)
abs_wf_t = np.abs(wf_t)
print("abs max: ", abs_wf_t.max())
print("abs min: ", abs_wf_t.min())

print("wf_t[0]: \n", wf_t[0])
print("wf_t[1]: \n", wf_t[1])

