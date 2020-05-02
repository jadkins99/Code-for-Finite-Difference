import numpy as np
def u_p_inf_call(xgrid,tau,k):
     # page 46 and 166 plus some algebra
    return np.exp(((1/2)*(k+1)*xgrid) - (1/4)*((k+1)**2)*tau)
    #return np.exp(0.25*((k+1)**2)*tau)*np.maximum(np.exp(0.5*(k+1)*xgrid-np.exp(0.5*(k-1)*xgrid)),0)