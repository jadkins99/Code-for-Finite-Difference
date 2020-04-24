from explicit_fd import explicit_fd
from call_payoff import call_payoff
from put_payoff import  put_payoff
from u_m_inf_call import u_m_inf_call
from u_p_inf_call import u_p_inf_call
from u_m_inf_put import u_m_inf_put
from u_p_inf_put import u_p_inf_put
import  numpy as np
import matplotlib.pyplot as plt

# Values chosen to duplicate figure 8.6 on page 143

E = 10
r = 0.05
sigma = 0.2
S_max = 16.0
T = 0.5 # six months to expiry
a = 0.25
Nx = 1000
k = r/(0.5*sigma**2) # page 136
''' change variables from  financial variables to dimensionless
 x = ln(S/E)
 tau = 0.5*(sigma**2)*(T-t)
v = P/E page 76 and 122
U = v*exp(-alpha*x-beta*tau) '''





uMatrix,xgrid = explicit_fd(E,r,sigma,T,S_max,Nx,a,put_payoff,u_m_inf_put,u_p_inf_put)

# transform variables into financial variables
S = E*np.exp(xgrid)

t = T

tau = 0.5*(sigma**2)*(T-t)

# page 136

#sPow = S**(0.5*(1-k))
#sMat = np.tile(sPow,(int(M),1))

V = np.exp(-1*0.5*(k-1)*xgrid - 0.25*((k+1)**2)*tau)*uMatrix
#V = (E**(0.5*(1+k)))*(xgrid**(0.5*(1-k)))*(np.exp(0.125*((k+1)**2)*(sigma**2)*(T-t)))*uMatrix
C = E*V
plt.plot(S,C[-1,:])
plt.show()



