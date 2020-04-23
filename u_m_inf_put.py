import numpy as np
def u_m_inf_put(x,tau,k):
    
    return np.exp(0.25*((k+1)**2)*tau)*max(np.exp(0.25*((k-1)**2)*x) - np.exp(0.25*((k+1)**2)*x),0.0)
