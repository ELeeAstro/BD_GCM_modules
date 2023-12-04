import numpy as np

qc = 1e-3
Nc = 5e8
rhoc = 3950.0
sig = 1

expterm = np.exp(9.0/2.0 * sig**2)

a3 = (3.0*qc)/(4.0*np.pi*rhoc*Nc*expterm)
a = np.cbrt(a3)

print(a*1e6)