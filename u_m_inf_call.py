import numpy as np
def u_m_inf_call(xgrid,tau,k):
    # page 46 and 166
    return 0.0
    #return np.exp(0.25*((k+1)**2)*tau)*np.maximum(np.exp(0.5*(k+1)*xgrid-np.exp(0.5*(k-1)*xgrid)),0)