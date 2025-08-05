import numpy as np
from main import *

## Expected β = 4

def gfunction(x, d):
    x1 = x[0]
    x2 = x[1]
    
    term1 = 8 * np.exp(-((x1 + 1)**2 + (x2 + 1)**2))
    term2 = 2 * np.exp(-((x1 - 5)**2 + (x2 - 4)**2))
    term3 = (x1 * x2) / 10

    g = 7 - term1 + term2 + 1 + term3
    return g

xvar = [
    {'varname': 'x1', 'vardist': 'normal', 'varmean': 2, 'varstd': 1},
    {'varname': 'x2', 'vardist': 'normal', 'varmean': 2, 'varstd': 1},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

example = Reliability(xvar, dvar, gfunction)
example.sampling_enhanced(100000, 5000, 0.1)