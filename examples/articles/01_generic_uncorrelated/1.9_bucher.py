from struct_reliability import *

## Expected β = 5.11

def gfunction(x, d):
    x1 = x[0]
    x2 = x[1]
    g = x1 * x2 - 146.14
    return g

xvar = [
    {'varname': 'x1', 'vardist': 'normal', 'varmean': 7.80644e4, 'varstd': 1.17097e4},
    {'varname': 'x2', 'vardist': 'normal', 'varmean': 0.0104, 'varstd': 0.00156},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

example = Reliability(xvar, dvar, gfunction)
example.bucher(1000, 5000, 0.01)
