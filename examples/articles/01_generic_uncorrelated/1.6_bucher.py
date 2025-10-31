from struct_reliability import *

## Expected β = 4.0

def gfunction(x, d):  
  g = (120/x[0]) - (x[1]/x[0]) - 1
  return g

xvar = [
  {'varname': 'x1', 'vardist': 'normal', 'varmean': 4, 'varstd': 1},
  {'varname': 'x2', 'vardist': 'normal', 'varmean': 4, 'varstd': 1},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

example = Reliability(xvar, dvar, gfunction)
example.bucher(1000, 5000, 0.01)