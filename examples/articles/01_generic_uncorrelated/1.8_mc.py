from struct_reliability import *

## Expected β = 2.90

def gfunction(x, d):
    x1 = x[0]
    x2 = x[1]
    g = x1**4 + 2 * x2**4 - 20
    return g

xvar = [
    {'varname': 'x1', 'vardist': 'normal', 'varmean': 10, 'varstd': 5},
    {'varname': 'x2', 'vardist': 'normal', 'varmean': 10, 'varstd': 5},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

example = Reliability(xvar, dvar, gfunction)
example.mc(1000, 5000, 0.01)
