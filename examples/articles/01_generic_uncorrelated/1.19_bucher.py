from struct_reliability import *
from scipy.stats import norm
import numpy as np

## Expected β ≈ 3.73 (≈ -> less reliable answer)

def gfunction(x, d):
    x1 = x[0]
    x2 = x[1]
    g = -0.16 * (x1 - 1)**3 - x2 + 4 - 0.04 * np.cos(x1 * x2 )
    return g

xvar = [
    {'varname': 'x1', 'vardist': 'normal', 'varmean': 0.0, 'varstd': 1.0},
    {'varname': 'x2', 'vardist': 'normal', 'varmean': 0.0, 'varstd': 1.0},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

example = Reliability(xvar, dvar, gfunction)
example.bucher(1000, 5000, 0.01)
