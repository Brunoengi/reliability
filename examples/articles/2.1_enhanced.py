import numpy as np
from main import *

## Expected β = 1.85

def gfunction(x, d):
  g = 3 * x[3] - (2 * x[4]) / (x[1] + x[2]) * np.sin(x[5] / 2 * np.sqrt((x[1] + x[2]) / x[0]))
  return g

xvar = [
  {'varname': 'x1', 'vardist': 'lognormal', 'varmean': 1, 'varstd': 0.05},
  {'varname': 'x2', 'vardist': 'lognormal', 'varmean': 1, 'varstd': 0.1},
  {'varname': 'x3', 'vardist': 'lognormal', 'varmean': 0.1, 'varstd': 0.01},
  {'varname': 'x4', 'vardist': 'lognormal', 'varmean': 0.5, 'varstd': 0.05},
  {'varname': 'x5', 'vardist': 'lognormal', 'varmean': 1, 'varstd': 0.2},
  {'varname': 'x6', 'vardist': 'lognormal', 'varmean': 1, 'varstd': 0.2},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00}
]

example = Reliability(xvar, dvar, gfunction)
example.mc(100, 5000, 0.01)