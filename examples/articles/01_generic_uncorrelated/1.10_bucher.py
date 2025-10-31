from struct_reliability import *

## Expected β = 5.21

def gfunction(x, d):
    x1 = x[0]
    x2 = x[1]
    g = x1 * x2 - 1140
    return g

xvar = [
    {'varname': 'x1', 'vardist': 'lognormal', 'varmean': 38.0, 'varstd': 3.8},
    {'varname': 'x2', 'vardist': 'lognormal', 'varmean': 54.0, 'varstd': 2.7},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

example = Reliability(xvar, dvar, gfunction)
example.bucher(1000, 5000, 0.01)
