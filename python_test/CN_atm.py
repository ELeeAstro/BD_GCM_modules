import numpy as np
import matplotlib.pylab as plt
import time

def  solve_tridiag(a,b,c,d):
     
#	 a - sub-diagonal (means it is the diagonal below the main diagonal)
#	 b - the main diagonal
#	 c - sup-diagonal (means it is the diagonal above the main diagonal)
#	 d - right part
#	 x - the answer

# initialize c-prime and d-prime
  n = len(d)
  cp = np.zeros(n)
  dp = np.zeros(n)

  cp[0] = c[0]/b[0]
  dp[0] = d[0]/b[0]
# solve for vectors c-prime and d-prime
  for i in range(1,n):
    m = b[i]-cp[i-1]*a[i]
    cp[i] = c[i]/m
    dp[i] = (d[i]-dp[i-1]*a[i])/m
 
  x = np.zeros(n)
# initialize x
  x[-1] = dp[-1]
#solve for x from the vectors c-prime and d-prime
  for i in range(n-2,-1,-1):
    x[i] = dp[i]-cp[i]*x[i+1]
  return x

kb = 1.380649e-23

nlay = 60
nlev = nlay + 1

Rd_air = 3571.0 
grav = 316.0

T = np.zeros(nlay)
T[:] = 1000.0

pe = np.logspace(-4,3,nlev) * 1e5
pl = np.zeros(nlay)
pl[:] = (pe[1:] - pe[:-1]) / np.log(pe[1:]/pe[:-1])

Kzz_e = np.zeros(nlay-1)
Kzz_e[:] = 1e4

nd_e = np.zeros(nlay)
nd_e[:] = pe[1:]/(kb * T[0]) 

alte = np.zeros(nlev)
delz = np.zeros(nlay)
alte[-1] = 0.0
for k in range(nlev-2,-1,-1):
  alte[k] = alte[k+1] + (Rd_air*T[k])/grav * np.log(pe[k+1]/pe[k])
  delz[k] = alte[k] - alte[k+1]

delz_mid = np.zeros(nlay-1)
for k in range(nlay-1): 
  delz_mid[k] = (alte[k] + alte[k+1])/2.0 - (alte[k+1] + alte[k+2])/2.0
  delz_mid[-1] = delz[-1]

print(delz[:])
print(delz_mid[:])


U = np.zeros(nlay)
U[:] = 1e-30
U[-1] = 1e-5

t_end = 1000000.0
t_now = 0
dt = 1000.0

a = np.zeros(nlay-2)
b = np.zeros(nlay-1)
c = np.zeros(nlay-2)


alp = np.zeros(nlay-1)
alp[:] = (Kzz_e[:]*dt)/(delz_mid[:]**2)

a[:] = -1.0
b[:] = 2.0*(1.0/alp[:] + 1.0)
c[:] = -1.0

d = np.zeros(nlay)

figure, ax = plt.subplots()

plot1, = ax.plot(U[:], pl[:]/1e5)

plt.yscale('log')
plt.xscale('log')

plt.xlabel('q',fontsize=16)
plt.ylabel('p [bar]',fontsize=16)

plt.gca().invert_yaxis()

plt.ion()

while (t_now < t_end):

  #d[0] = 2.0*alp*U[0] + 2.0*(1.0 - alp)*U[1] + alp*U[2]
  for k in range(1,nlay-1):
    d[k] = U[k-1] + 2.0*(1.0/alp[k] - 1.0)*U[k] + U[k+1]
  #d[-1] = alp*U[-3] + 2.0*(1.0 - alp)*U[-2] + 2.0*alp*U[-1]

  #print('d',d[:])

  U[1:-1] = solve_tridiag(a,b,c,d[1:-1])

  U[0] = U[1]

  #print('u',U[:])

  plot1.set_xdata(U[:])
  plot1.set_ydata(pl[:]/1e5)

  plt.title(str(t_now))
    
  figure.canvas.draw()
  figure.canvas.flush_events()
  time.sleep(0.01)

  plt.show()

  #U[-1] = 1.0
  #U[0] = U[1]

  t_now = t_now + dt

