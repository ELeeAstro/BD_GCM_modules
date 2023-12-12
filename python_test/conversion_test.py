import numpy as np


amu = 1.66053906660e-24
a = 1.0e-2 * 1e-4

bulk_den = 4.23
mol_wght = 79.866
mass = mol_wght * amu
dV = mass / bulk_den
Nl = (4.0/3.0 * np.pi * a**3) / dV

print(Nl)