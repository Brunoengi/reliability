from main import *

## Expected β = 2.53

def gfunction(x, d):  
  g = x[0]**3 + x[1]**3 - 18
  return g

xvar = [
  {'varname': 'x1', 'vardist': 'normal', 'varmean': 10, 'varstd': 5},
  {'varname': 'x2', 'vardist': 'normal', 'varmean': 9.9, 'varstd': 5},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

example = Reliability(xvar, dvar, gfunction)
example.mc(100000, 5000, 0.01)