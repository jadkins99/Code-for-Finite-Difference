import numpy as np
def call_payoff(xgrid,k):
    
    nx = len(xgrid)
    # page 77
    return np.maximum(np.exp(0.5*(k+1)*xgrid)-np.exp(0.5*(k-1)*xgrid),np.zeros(nx))
    
    