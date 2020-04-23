import numpy as np
def u_p_inf_call(x,tau,k):
    
    return np.exp(0.25*((k+1)**2)*tau)*max(np.exp(0.25*((k+1)**2)*x) - np.exp(0.25*((k-1)**2)*x),0.0)
