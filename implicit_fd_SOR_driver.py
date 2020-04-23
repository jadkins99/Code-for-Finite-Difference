import numpy as np
import matplotlib.pyplot as plt
from implicit_fd_SOR import implicit_fd_SOR

E = 10.0
r = 0.05
sigma = 0.2
T = 0.5 # six months to expiry
k = r/(0.5*sigma**2) # page 136
t = 0

''' change variables from  financial variables to dimensionless
 x = ln(S/E)
 tau = 0.5*(sigma**2)*(T-t)
v = P/E page 76 and 122
U = v*exp(-alpha*x-beta*tau) '''
S_min = 0.00000001
S_max = 160.0

x_min = np.log(S_min/E)
x_max = np.log(S_max/E)

nx = 1000

dx = (x_max-x_min)/nx
a = 0.25

dt = a*dx**2 # page 140

tau_max = 0.5*(sigma**2)*T
M = np.ceil(tau_max/dt) # page 141


u,xgrid = implicit_fd_SOR(r,sigma,x_min,x_max,nx,tau_max,M)

S = E*np.exp(xgrid)
tau = 0.5*(sigma**2)*(T-t)

sPow = S**(0.5*(1-k))
sMat = np.tile(sPow,(int(M),1))

V    = (E**(0.5*(1+k))) * sMat * np.exp((-1*(1/4)*((k+1)**2)*tau))*u

plt.plot(S,V[-1,:])
plt.show()
