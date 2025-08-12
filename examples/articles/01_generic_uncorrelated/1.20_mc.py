from main import *

## Expected β = 3.43

def gfunction(x, d):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    g = (1 / 40) * x1**4 + 2 * x2**2 + x3 + 3
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