import numpy as np
from struct_reliability import *

## Expected β = 3.633

def gfunction(x, d):
    g = 2 - np.exp(x[4] * x[2] / x[0]) + ((np.exp(x[4])- 2) / (np.exp(-x[5]) - 1)) * (np.exp(-x[5] * x[2] / x[0]) - 1) - x[3] / x[1]
    return g

xvar = [
  {'varname': 'x1', 'vardist': 'lognormal', 'varmean': 5490, 'varstd': 1098},
  {'varname': 'x2', 'vardist': 'lognormal', 'varmean': 17100, 'varstd': 3420},
  {'varname': 'x3', 'vardist': 'lognormal', 'varmean': 549, 'varstd': 109.8},
  {'varname': 'x4', 'vardist': 'lognormal', 'varmean': 4000, 'varstd': 800},
  {'varname': 'x5', 'vardist': 'normal', 'varmean': 0.42, 'varstd': 0.084},
  {'varname': 'x6', 'vardist': 'normal', 'varmean': 6, 'varstd': 1.2},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00}
]

example = Reliability(xvar, dvar, gfunction)
example.sampling_enhanced(100000, 5000, 0.01)