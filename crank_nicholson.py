import numpy as np

def crank_nicholson(E,r,sigma,T,s_max,Nx,M,pay_off,u_m_inf,u_p_inf):
    
    s_min = 0.00000001

    #x = log(S/E)
    xLeft = np.log(s_min/E)
    xRight = np.log(s_max/E)

    # tau = 1/2sigma^2(T-t)
    tau_max = 0.5*(sigma**2)*T

    dx = (xRight-xLeft)/Nx
    dt = tau_max/M
    a = dt/(dx**2) # page 140
    a2 = a/2
    print("a = ",a)
    xgrid = np.linspace(xLeft,xRight,Nx)

    k = r/(0.5*sigma**2) # page 136
  
  
    # initial conditions
    tau = 0.0
    values = pay_off(xgrid,k)
    
    MSize = int(Nx-2)

    # eq. 8.30 and 8.31 on page 157
    Cmat = (1+a)*np.eye(MSize,MSize,k=0) + (-a2)*np.eye(MSize,MSize,k=1) + (-a2)*np.eye(MSize,MSize,k=-1) 

    for i in range(1,int(M)):
        tau = i*dt

        b = (1-a)*values[1:-1] + a2*(values[2:]+values[0:-2])
        values[0] = u_m_inf(xgrid[0],tau,k)
        values[-1] = u_p_inf(xgrid[-1],tau,k)

        b[0] += a2*values[0]
        b[-1] += a2*values[-1]
        # eq. 8.32 and eq. 8.34 on page 157    
        newu = np.linalg.solve(a = Cmat,b = b)
        values[1:-1] = newu

    return values,xgrid

