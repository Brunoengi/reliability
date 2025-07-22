import numpy as np
from main import *

## Expected β = 2.07

def gfunction(x, d=None):
    x1 = x[0]
    x2 = x[1]

    term1 = 2.2257
    term2 = (0.025 * np.sqrt(2) / 27) * (x1 + x2 - 20)**3
    term3 = (33 / 140) * (x1 - x2)

    g = term1 - term2 + term3
    return g

xvar = [
    {'varname': 'x1', 'vardist': 'normal', 'varmean': 10, 'varstd': 3},
    {'varname': 'x2', 'vardist': 'normal', 'varmean': 10, 'varstd': 3},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

example = Reliability(xvar, dvar, gfunction)
example.sampling_enhanced(100000, 5000, 0.01)