from explicit_fd import explicit_fd
import  numpy as np
import matplotlib.pyplot as plt

# Values chosen to duplicate figure 8.6 on page 143

E = 10
r = 0.05
sigma = 0.2
T = 0.5 # six months to expiry
k = r/(0.5*sigma**2) # page 136


''' change variables from  financial variables to dimensionless
 x = ln(S/E)
 tau = 0.5*(sigma**2)*(T-t)
v = P/E page 76 and 122
U = v*exp(-alpha*x-beta*tau) '''
S_min = 0.00000001
S_max = 16.0

xLeft = np.log(S_min/E)
XRight = np.log(S_max/E)

Nx = 2000
dx = (XRight-xLeft)/Nx
a = 0.25
dt = a*dx**2 # page 140

tau_max = 0.5*(sigma**2)*T
M = np.ceil(tau_max/dt) # page 141

u,xgrid = explicit_fd(r,sigma,xLeft,XRight,Nx,tau_max,M)

# transform variables into financial variables
S = E*np.exp(xgrid)

t = 0

tau = 0.5*(sigma**2)*(T-t)

# page 136

sPow = S**(0.5*(1-k))
sMat = np.tile(sPow,(int(M),1))

V    = (E**(0.5*(1+k))) * sMat * np.exp((-1*(1/4)*((k+1)**2)*tau))*u

plt.plot(S,V[-1,:])
plt.show()



