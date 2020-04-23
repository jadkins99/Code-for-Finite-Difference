from u_m_inf_call import u_m_inf_call
from u_p_inf_call import u_p_inf_call
from u_p_inf_put import u_p_inf_put
from u_m_inf_put import u_m_inf_put
from put_payoff_transformed import  put_payoff
from call_payoff_transformed import call_payoff
import numpy as np


def explicit_fd(r,sigma,xLeft,xRight,Nx,tau_max,M):

    
    xgrid = np.linspace(xLeft,xRight,Nx)
    dx = (xRight-xLeft)/Nx
    dt = tau_max/Nx
    a = dt/(dx**2)
    k = r/(0.5*sigma**2)
    
    # check for stability
    if a >= 0.5:
        print("a is to large. This method will become unstable.")
        return False
    
    # initial conditions
    tau = 0.0
    oldu = call_payoff(xgrid,tau,k)

   
    un = np.zeros((int(M),int(Nx)))
    un[0,:] = oldu
   
    newu = np.zeros((int(Nx)))
    
    for m in range(1,int(M)):
        
        # update endpoints 
        newu[0] = u_m_inf_call(xgrid[0],tau,k)
        newu[-1] = u_p_inf_call(xgrid[-1],tau,k)

        # update newu
        newu[2:-1] = oldu[2:-1] + a*(oldu[1:-2] - 2*oldu[2:-1] + oldu[3:])

        # prepare for new iteration
        oldu = newu

        un[m,:] = newu
    
    return un,xgrid

