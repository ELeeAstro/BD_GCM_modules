import numpy as np
import matplotlib.pylab as plt

fname = 'bforce.txt'


data = np.loadtxt(fname)

lon = data[0,1:]
lat = data[1:,0]

print(lon, lat)

bforce = data[1:,1:]

print(bforce)

fig = plt.figure()

plt.contourf(lon,lat,bforce)
cbar = plt.colorbar()
cbar.solids.set_edgecolor("face")

plt.show()