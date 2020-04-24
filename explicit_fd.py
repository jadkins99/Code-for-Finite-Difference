import numpy as np


def explicit_fd(E,r,sigma,T,s_max,Nx,a,pay_off,u_m_inf,u_p_inf):

    s_min = 0.00000001

    #x = log(S/E)
    xLeft = np.log(s_min/E)
    xRight = np.log(s_max/E)

    # tau = 1/2sigma^2(T-t)
    tau_max = 0.5*(sigma**2)*T

    dx = (xRight-xLeft)/Nx
    dt = a*dx**2 # page 140

    xgrid = np.linspace(xLeft,xRight,Nx)

    dx = (xRight-xLeft)/Nx
    dt = tau_max/Nx
    a = dt/(dx**2)
    k = r/(0.5*sigma**2) # page 136
    M = np.ceil(tau_max/dt) # page 141  
  
    # initial conditions
    tau = 0.0
    oldu = pay_off(xgrid,k)

   
    uMat = np.zeros((int(M),int(Nx)))
    uMat[0,:] = oldu.copy()
   
    newu = np.zeros((int(Nx)))
    
    for m in range(1,int(M)):

        tau = m*dt
        
        # update endpoints 
        newu[0] = u_m_inf(xgrid[0],tau,k)
        newu[-1] = u_p_inf(xgrid[-1],tau,k)

        # update newu
        for n in range(1,len(newu)-1):
            newu[n] = oldu[n] + a*(oldu[n-1] - 2*oldu[n] + oldu[n+1])

        # prepare for new iteration
        oldu = newu.copy()

        uMat[m,:] = newu.copy()
    
    return uMat,xgrid

