import numpy as np

def implicit_fd(E,r,sigma,T,s_max,Nx,M,pay_off,u_m_inf,u_p_inf):

    s_min = 0.00000001

    #x = log(S/E)
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

    
    uMat = np.zeros((int(M),int(Nx)))
    uMat[0,:] = oldu

    newu = np.zeros((int(Nx)))

    b = np.zeros((int(Nx)))
    MSize = int(Nx)
    Mmat = (1+2*a)*np.eye(MSize,MSize,k=0) + (-a)*np.eye(MSize,MSize,k=1) + (-a)*np.eye(MSize,MSize,k=-1) 
    for i in range(1,int(M)):
        tau = i*dt
        
        b[0:] = oldu[0:].copy()
        oldu[0] = u_m_inf(xgrid[0],tau,k)
        oldu[-1] = u_m_inf(xgrid[-1],tau,k)

        b[0] += a*oldu[0]
        b[-1] += a*oldu[-1]
       
        newu = np.linalg.solve(a = Mmat,b = b)
        
        oldu = newu.copy()
        uMat[i,:] = newu.copy()

    return uMat,xgrid