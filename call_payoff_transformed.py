import numpy as np
def call_payoff(xgrid,tau,k):
    nx = len(xgrid)
   
    g = np.exp(0.25*((k+1)**2)*tau)*np.maximum(np.exp(0.5*(k+1)*xgrid)-np.exp(0.5*(k-1)*xgrid),np.zeros(nx))
    print(g)
    
    return g