import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

nlay = 60         
nsp = 13
ncld = 1

data1 = np.loadtxt('FMS_RC_ic.out',skiprows=1)
Pic = data1[:,1]
Tic = data1[:,2]
muic = data1[:,3]
qic = data1[:,4:]

data2 = np.loadtxt('FMS_RC_pp_2.out',skiprows=1)
Ppp = data2[:,1]
Tpp = data2[:,2]
rdbar = data2[:,3]
cpbar = data2[:,4]
kprime = data2[:,5] 
Kzz = data2[:,6] 


fig = plt.figure()

plt.plot(rdbar,Ppp/1e5,ls='dashed',lw=3,label='rd bar',c='orange')

plt.ylabel('Pressure [bar]')
plt.xlabel('Rd bar')
plt.legend()

plt.yscale('log')
plt.gca().invert_yaxis()

fig = plt.figure()

plt.plot(cpbar,Ppp/1e5,ls='dashed',lw=3,label='cp bar',c='orange')

plt.ylabel('Pressure [bar]')
plt.xlabel('cp bar')
plt.legend()

plt.yscale('log')
plt.gca().invert_yaxis()

fig = plt.figure()

plt.plot(kprime,Ppp/1e5,ls='dashed',lw=3,label='k prime',c='orange')

plt.ylabel('Pressure [bar]')
plt.xlabel('k prime')
plt.legend()

plt.yscale('log')
plt.gca().invert_yaxis()

fig = plt.figure()

plt.plot(Kzz*1e4,Ppp/1e5,ls='dashed',lw=3,label='Kzz',c='orange')

plt.ylabel('Pressure [bar]')
plt.xlabel('Kzz')
plt.legend()

plt.yscale('log')
plt.xscale('log')
plt.gca().invert_yaxis()

plt.xlim(1e1,1e12)

plt.show()
