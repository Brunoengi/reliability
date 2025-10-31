from struct_reliability import *
import numpy as np
from scipy.special import gamma

## Expected β = 1.734

def gfunction(x, d):
    g = x[0] - x[1]
    return g

xvar = [
  {'varname': 'x1', 'vardist': 'uniform', 'parameter1': 9, 'parameter2': 11},
  {'varname': 'x2', 'vardist': 'normal', 'varmean': 8, 'varstd': 1},

]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00}
]

example = Reliability(xvar, dvar, gfunction)
example.mc(1000, 5000, 0.01)