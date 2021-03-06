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
Nx = 100
M = 100
k = r/(0.5*sigma**2) # page 136


u,xgrid = explicit_fd(E,r,sigma,T,S_max,Nx,M,call_payoff,u_m_inf_call,u_p_inf_call)
#u,xgrid = implicit_fd(E,r,sigma,T,S_max,Nx,M,call_payoff,u_m_inf_call,u_p_inf_call)
#u,xgrid = crank_nicholson(E,r,sigma,T,S_max,Nx,M,call_payoff,u_m_inf_call,u_p_inf_call)

# transform variables into financial variables
S = E*np.exp(xgrid)

t = 0

# page 136

V = (E**(0.5*(1+k)))*(S**(0.5*(1-k)))*(np.exp(0.125*((k+1)**2)*(sigma**2)*(T-t)))*u

plt.plot(S,V)
plt.xlabel("S")
plt.ylabel("V")
plt.show()





