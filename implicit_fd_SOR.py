from u_m_inf_call import u_m_inf_call
from u_p_inf_call import u_p_inf_call
from u_p_inf_put import u_p_inf_put
from u_m_inf_put import u_m_inf_put
from put_payoff_transformed import  put_payoff
from call_payoff_transformed import call_payoff
import numpy as np
from SOR_solver import SOR_solver

def implicit_fd_SOR(r,sigma,xLeft,xRight,Nx,tau_max,M):

    xgrid = np.linspace(xLeft,xRight,Nx)
    dx = (xRight-xLeft)/Nx
    dt = tau_max/Nx
    a = dt/(dx**2)
    k = r/(0.5*sigma**2)
    omega = 1.0
    domega = 0.05
    old_loops = 10000
    eps = 0.00000001
    tau = 0

    oldu = call_payoff(xgrid,tau,k)
   
    u = np.zeros((int(M),int(Nx)))
    u[0,:] = oldu

    newu = np.zeros((int(Nx)))

    b = np.zeros((int(Nx)))

    for m in range(1,int(M)):

        tau = m*dt
        b[2:-1] = oldu[2:-1]
        
        # update endpoints 
        oldu[0] = u_m_inf_call(xgrid[0],tau,k)
        oldu[-1] = u_p_inf_call(xgrid[-1],tau,k)
        
        new_u,loops = SOR_solver(oldu,b,a,omega,eps)
      
        if loops>old_loops:
            print("HELLLOOOO")            
            domega *= -1

        omega+=domega

        old_loops = loops

        oldu = new_u

        u[m,:] = newu

    return u,xgrid