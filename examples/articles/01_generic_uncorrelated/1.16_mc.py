from struct_reliability import *

## Expected β = 1.819

def gfunction(x, d):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    g = 2 - x2 - 0.1 * x1**2 + 0.06 * x3**3
    return g

xvar = [
    {'varname': 'x1', 'vardist': 'normal', 'varmean': 0.0, 'varstd': 1.0},
    {'varname': 'x2', 'vardist': 'normal', 'varmean': 0.0, 'varstd': 1.0},
    {'varname': 'x3', 'vardist': 'normal', 'varmean': 0.0, 'varstd': 1.0},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

example = Reliability(xvar, dvar, gfunction)
example.mc(1000, 5000, 0.01)