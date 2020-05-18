import numpy as np
def u_m_inf_put(xgrid,tau,k):
    # page 46 and 166 plus some algebra
    return np.exp(((1/2)*(k-1)*xgrid) - (1/4)*((k+1)**2)*tau)
    