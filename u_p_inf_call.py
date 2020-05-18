import numpy as np
def u_p_inf_call(xgrid,tau,k):
     # page 46 and 166 plus some algebra

    return np.exp(-0.25*((k+1)**2)*tau)*(np.exp(0.5*xgrid*(k+1))-np.exp(0.5*xgrid*(k-1)))    