import numpy as np
import scipy.sparse as sp

def crank_nicholson(E,r,sigma,T,s_max,Nx,a,pay_off,u_m_inf,u_p_inf):


    
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
    uMat[0,:] = oldu

    newu = np.zeros((int(Nx)))

    b = np.zeros((int(Nx)))
    MSize = int(Nx)

    a2 = a/2
    Cmat = (1+a)*np.eye(MSize,MSize,k=0) + (-a2)*np.eye(MSize,MSize,k=1) + (-a2)*np.eye(MSize,MSize,k=-1) 



    for i in range(1,int(M)):
        tau = i*dt
        
        for j in range(1,Nx-1):
            b[j] = (1-a)*oldu[j]+a2*(oldu[j+1]+oldu[j-1])

        

        oldu[0] = u_m_inf(xgrid[0],tau,k)
        oldu[-1] = u_m_inf(xgrid[-1],tau,k)

        b[1] += a2*oldu[0]
        b[-2] += a2*oldu[-1]
        
        newu = np.linalg.solve(a = Cmat,b = b)
        
       # print(len(newu))
        oldu = newu.copy()
        uMat[i,:] = newu.copy()

    return uMat,xgrid

