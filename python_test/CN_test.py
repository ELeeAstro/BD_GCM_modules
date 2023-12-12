import numpy as np
import matplotlib.pylab as plt

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

rho = 2.7
C = 0.2174
k = 0.49
alpha = k/(rho*C)

Ll = 0
Lu = 10
dt = 0.2
dx = 1.0

nlay = Lu

alp = (alpha*dt)/(dx**2)

U = np.zeros(nlay)
U[0] = 100
U[-1] = 100

t_end = 10.0
t_now = 0
dt = 0.2

a = np.zeros(nlay)
b = np.zeros(nlay)
c = np.zeros(nlay)

a[:] = -1.0
b[:] = 2.0*(1.0/alp + 1.0)
c[:] = -1.0


d = np.zeros(nlay)

for t in range(1111):
  if (t_now > t_end):
    break

  #d[0] = 2.0*alp*U[0] + 2.0*(1.0 - alp)*U[1] + alp*U[2]
  for k in range(1,nlay-1):
    d[k] = U[k-1] + 2.0*(1.0/alp - 1.0)*U[k] + U[k+1]
  #d[-1] = alp*U[-3] + 2.0*(1.0 - alp)*U[-2] + 2.0*alp*U[-1]

  print('d',d[:])

  U[1:-1] = solve_tridiag(a,b,c,d[1:-1])

  print('u',U[:])

  #U[-1] = 1.0
  #U[0] = U[1]

  t_now = t_now + dt