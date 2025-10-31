from struct_reliability import *

## Expected β = 2.41

def gfunction(x, d):
  g = 18.461 - 7.477 * (10**10 * (x[0]/(x[1]**3)))

  return g

xvar = [
  {'varname': 'x1', 'vardist': 'normal', 'varmean': 0.001, 'varstd': 0.0002},
  {'varname': 'x2', 'vardist': 'normal', 'varmean': 250, 'varstd': 37.5},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00}
]

example = Reliability(xvar, dvar, gfunction)
example.sampling_enhanced(1000, 5000, 0.01)