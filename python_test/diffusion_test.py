import numpy as np
import matplotlib.pylab as plt


nt = 100

t_end = 60.0

nlay = 60
nlev = nlay + 1
plev = np.logspace(-6,3,nlev)
play = np.zeros(nlay)
for  k in range(nlay):
  play[k] = (plev[k+1] - plev[k]) / (np.log(plev[k+1]/plev[k]))
print(plev)
print(play)

plev = plev * 1e5
play = play * 1e5

Rd = 3714.0
grav = 1000.0

Tlay = np.zeros(nlay)
Tlay[:] = 1000.0

# Calculate vertical height grid
alt_lev = np.zeros(nlev)
alt_lev[-1] = 0.0
for k in range(nlev-2, -1, -1):
  alt_lev[k] = alt_lev[k+1] + (Rd*Tlay[k])/grav * np.log(plev[k+1]/plev[k])


alt_lev[:] = alt_lev[::-1]
print(alt_lev)

rho = np.zeros(nlay)
rho[:] = play[:]/(Rd * Tlay[:])

Kzz = np.zeros(nlay)
Kzz[:] = 1e7

Kzz_lev = np.zeros(nlay)
Kzz_lev[:] = 1e7
Kzz_lev[:] = Kzz_lev[::-1]

q = np.zeros(nlay)
q[-1] = 1.0
q[:] = q[::-1]

# Prepare grid
hr = np.zeros(nlay)
hl = np.zeros(nlay)
d1l = np.zeros(nlay)
d1m = np.zeros(nlay)
d1r = np.zeros(nlay)
d2l = np.zeros(nlay)
d2m = np.zeros(nlay)
d2r = np.zeros(nlay)

#--- compute grid point differences ---
for i in range(nlay):
  hr[i] = alt_lev[i+1] - alt_lev[i]
#print(hr[:])
for i in range(1,nlay):
  hl[i] = alt_lev[i] - alt_lev[i-1]
#print(hl[:])
for i in range(1,nlay):
  d1l[i] = -hr[i]/((hr[i]+hl[i])*hl[i])
  d1m[i] =  (hr[i]-hl[i])/(hl[i]*hr[i])
  d1r[i] =  hl[i]/((hr[i]+hl[i])*hr[i])
  d2l[i] =  2.0/((hr[i]+hl[i])*hl[i])
  d2m[i] = -2.0/(hr[i]*hl[i])
  d2r[i] =  2.0/((hr[i]+hl[i])*hr[i])

h1 = alt_lev[1]-alt_lev[0]
h2 = alt_lev[2]-alt_lev[0]
d1l[0] = -(h1+h2)/(h1*h2)
d1m[0] =  h2/(h1*(h2-h1))
d1r[0] = -h1/(h2*(h2-h1))

h1 = alt_lev[-2]-alt_lev[-1]
h2 = alt_lev[-3]-alt_lev[-1]
d1r[-1] = -(h1+h2)/(h1*h2)
d1m[-1] =  h2/(h1*(h2-h1))
d1l[-1] = -h1/(h2*(h2-h1))

rate = np.zeros(nlay)

for t in range(nt):
  dt = 9.0e99
  for i in range(1,nlay):
    dt = np.minimum(dt,0.33*(alt_lev[i]-alt_lev[i-1])**2/Kzz[i])

  t_now = 0.0
  for k in range(1000000):

    if (t_now >= t_end):
      break

    influx = -Kzz_lev[0]*(d1l[0]*q[0] + d1m[0]*q[1] + d1r[0]*q[2])
    q[-1] = (-d1l[-1]*q[-3] - d1m[-1]*q[-2])/d1r[-1]

    for i in range(1,nlay-1):
      d1 = d1l[i]*q[i-1] + d1m[i]*q[i] + d1r[i]*q[i+1]
      d2 = d2l[i]*q[i-1] + d2m[i]*q[i] + d2r[i]*q[i+1]

      d1nD = d1l[i]*Kzz_lev[i-1] + d1m[i]*Kzz_lev[i] + d1r[i]*Kzz_lev[i+1]

      rate[i] = Kzz_lev[i]*d2 + d1nD*d1
      #print(i, rate[i])

    for i in range(nlay):
      q[i] = q[i] + rate[i]*dt

    influx = -Kzz_lev[0]*(d1l[0]*q[0] + d1m[0]*q[1] + d1r[0]*q[2])
    q[-1] = (-d1l[-1]*q[-3] - d1m[-1]*q[-2])/d1r[-1]

    #print(t, t_now, dt)
    #print(q) 

    t_now = t_now + dt

q = q[::-1]

fig = plt.figure()

plt.plot(q,play/1e5)

plt.xscale('log')
plt.yscale('log')

plt.gca().invert_yaxis()

plt.show()
