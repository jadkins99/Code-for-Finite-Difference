from crank_nicholson import crank_nicholson
from call_payoff import call_payoff
from put_payoff import  put_payoff
from u_m_inf_call import u_m_inf_call
from u_p_inf_call import u_p_inf_call
from u_m_inf_put import u_m_inf_put
from u_p_inf_put import u_p_inf_put
import  numpy as np
import matplotlib.pyplot as plt

# Values chosen to duplicate figure 8.6 on page 143

E = 100
r = 0.05
sigma = 0.2
S_max = 200
T = 0.5 # six months to expiry
a = 0.25
Nx = 200
k = r/(0.5*sigma**2) # page 136
''' change variables from  financial variables to dimensionless
 x = ln(S/E)
 tau = 0.5*(sigma**2)*(T-t)
v = P/E page 76 and 122
U = v*exp(-alpha*x-beta*tau) '''





uMatrix,xgrid = crank_nicholson(E,r,sigma,T,S_max,Nx,a,call_payoff,u_m_inf_call,u_p_inf_call)

# transform variables into financial variables
S = E*np.exp(xgrid)

tau_max = 0.5*(sigma**2)*T
dt = tau_max/Nx
t = 0

tau = 0.5*(sigma**2)*(T-t)
M = np.ceil(tau_max/dt)
m = int(np.floor(tau/dt))

# page 136



#V = E*np.exp(-1*0.5*(k-1)*xgrid - 0.25*((k+1)**2)*tau)*uMatrix

V = (E**(0.5*(1+k)))*(S**(0.5*(1-k)))*(np.exp(0.125*((k+1)**2)*(sigma**2)*(T-t)))*uMatrix

if m == Nx:
    m -=1


#C = E*V
plt.plot(S,V[0,:])
plt.xlabel("S")
plt.ylabel("V")
plt.show()





