from crank_nicholson import crank_nicholson
from explicit_fd import explicit_fd
from implicit_fd import implicit_fd
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
S_max = 16
T = 0.5 # six months to expiry
Nx = 160
M = 1000
k = r/(0.5*sigma**2) # page 136


#uMatrix,xgrid = crank_nicholson(E,r,sigma,T,S_max,Nx,M,put_payoff,u_m_inf_put,u_p_inf_put)
#uMatrix,xgrid = implicit_fd(E,r,sigma,T,S_max,Nx,M,call_payoff,u_m_inf_call,u_p_inf_call)
#uMatrix,xgrid = explicit_fd(E,r,sigma,T,S_max,Nx,M,call_payoff,u_m_inf_call,u_p_inf_call)
u,xgrid = explicit_fd(E,r,sigma,T,S_max,Nx,M,call_payoff,u_m_inf_call,u_p_inf_call)
# transform variables into financial variables
S = E*np.exp(xgrid)

#tau_max = 0.5*(sigma**2)*T
#dt = tau_max/M
t = 0

#tau = 0.5*(sigma**2)*(T-t)

#m = int(np.floor(tau/dt))

# page 136

V = (E**(0.5*(1+k)))*(S**(0.5*(1-k)))*(np.exp(0.125*((k+1)**2)*(sigma**2)*(T-t)))*u

#if m == M:
 #   m -=1

#plt.plot(S,V[-1,:])
plt.plot(S,V)
plt.xlabel("S")
plt.ylabel("V")
plt.show()





