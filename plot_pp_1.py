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

data2 = np.loadtxt('FMS_RC_pp_1.out',skiprows=1)
Ppp = data2[:,1]
Tpp = data2[:,2]
dTrad = data1[:,3]
dTconv = data1[:,4]
mupp = data2[:,5]
qpp = data2[:,6:] 


fig = plt.figure()

plt.plot(Tic,Pic/1e5,ls='dashed',lw=3,label='Initial Conditions',c='orange')
plt.plot(Tpp,Ppp/1e5,lw=3,label='Numerical Result',c='blue')

plt.ylabel('Pressure [bar]')
plt.xlabel('Temperature [K]')
plt.legend()

plt.yscale('log')
plt.gca().invert_yaxis()


fig = plt.figure()

plt.plot(dTrad,Ppp/1e5,ls='dashed',lw=3,label='dT rad',c='orange')
plt.plot(dTconv,Ppp/1e5,lw=3,label='dT conv',c='blue')

plt.ylabel('Pressure [bar]')
plt.xlabel('dT [K s-1]')
plt.legend()

plt.yscale('log')
plt.xscale('log')
plt.gca().invert_yaxis()

fig = plt.figure()

plt.plot(muic,Pic/1e5,ls='dashed',lw=3,label='mu ic',c='orange')
plt.plot(mupp,Ppp/1e5,lw=3,label='mu pp',c='blue')

plt.ylabel('Pressure [bar]')
plt.xlabel('mu [g mol$^{-1}$]')
plt.legend()

plt.yscale('log')
plt.gca().invert_yaxis()


fig = plt.figure()

lab = ['OH','H2','H2O','H','CO','CO2','O','CH4','C2H2','NH3','N2','HCN','He']

col = sns.color_palette("husl", nsp)

for i in range(nsp):
  plt.plot(qpp[:,i],Ppp/1e5,lw=2,label=lab[i],ls='solid',c=col[i])
for i in range(nsp):
  plt.plot(qic[:,i],Pic/1e5,lw=2,ls='dashed',c=col[i])

plt.ylabel('Pressure [bar]')
plt.xlabel('VMR')
plt.legend()

plt.yscale('log')
plt.xscale('log')
plt.gca().invert_yaxis()

fig = plt.figure()

plt.plot(qpp[:,nsp],Ppp/1e5,lw=2,ls='solid',c=col[0],label='qv')
plt.plot(qpp[:,nsp+1],Ppp/1e5,lw=2,ls='solid',c=col[1],label='qc')

plt.ylabel('Pressure [bar]')
plt.xlabel('q')
plt.legend()

plt.yscale('log')
plt.xscale('log')
plt.gca().invert_yaxis()

plt.show()
