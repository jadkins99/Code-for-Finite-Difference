import numpy as np
def put_payoff(xgrid,k):
    nx = len(xgrid)
    
    return np.maximum(np.exp(0.5*(k-1)*xgrid-np.exp(0.5*(k+1)*xgrid)),np.zeros(nx))