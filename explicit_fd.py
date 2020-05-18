import numpy as np

def explicit_fd(E,r,sigma,T,s_max,Nx,M,pay_off,u_m_inf,u_p_inf):

    s_min = 0.00000001

    #x = log_2(S/E)
    xLeft = np.log(s_min/E)
    xRight = np.log(s_max/E)

    # tau = 1/2sigma^2(T-t)
    tau_max = 0.5*(sigma**2)*T

    dx = (xRight-xLeft)/Nx
    dt = tau_max/M
    a = dt/(dx**2) # page 140
    print("a = ",a)
    xgrid = np.linspace(xLeft,xRight,Nx)
    
    
    k = r/(0.5*sigma**2) # page 136
    
  
    # initial conditions
    tau = 0.0
    oldu = pay_off(xgrid,k)

   
    newu = np.zeros((int(Nx)))
    values = np.zeros((int(Nx)))
    for m in range(1,int(M+1)):

        tau = m*dt
        
        # update endpoints 
        newu[0] = u_m_inf(xgrid[0],tau,k)
        newu[-1] = u_p_inf(xgrid[-1],tau,k)

        # update newu

        newu[1:-1] = oldu[1:-1] + a*(oldu[0:-2] - 2*oldu[1:-1] + oldu[2:])
        # prepare for new iteration
        oldu = newu.copy()

        values = oldu
        
    return values,xgrid

