import numpy as np
from struct_reliability import *

## Expected β = 2.05

def gfunction(x, d):
  #Sb
  Sb = x[5] / ((x[0] * (x[3] - 2 * x[2])**3) / (6 * x[3]) + (x[1] / (6 * x[3])) * (x[3]**3 - (x[3] - 2 * x[2])**3))
  
  #Ts
  Ts = x[5] / ((0.8 * x[1] * x[2]**2 + 0.4 * (x[0]**3 * (x[3] - 2 * x[2])) / x[2]))

  g = 460 - np.sqrt(Sb**2 + 3 * Ts**2)

  return g

xvar = [
  {'varname': 'x1', 'vardist': 'normal', 'varmean': 12, 'varstd': 0.06},
  {'varname': 'x2', 'vardist': 'normal', 'varmean': 65, 'varstd': 0.325},
  {'varname': 'x3', 'vardist': 'normal', 'varmean': 14, 'varstd': 0.07},
  {'varname': 'x4', 'vardist': 'normal', 'varmean': 85, 'varstd': 0.425},
  {'varname': 'x5', 'vardist': 'normal', 'varmean': 3.5e6, 'varstd': 1.75e5},
  {'varname': 'x6', 'vardist': 'normal', 'varmean': 3.1e6, 'varstd': 1.55e5},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00}
]

example = Reliability(xvar, dvar, gfunction)
example.sampling_enhanced(1000, 5000, 0.01)