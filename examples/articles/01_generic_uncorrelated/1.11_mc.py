from main import *

## Expected β = 2.184

def gfunction(x, d):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    g = x1 * x2 - 2000 * x3
    return g

xvar = [
    {'varname': 'x1', 'vardist': 'normal', 'varmean': 0.32, 'varstd': 0.032},
    {'varname': 'x2', 'vardist': 'normal', 'varmean': 1.4e6, 'varstd': 7.0e4},
    {'varname': 'x3', 'vardist': 'lognormal', 'varmean': 100.0, 'varstd': 40.0},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

example = Reliability(xvar, dvar, gfunction)
example.mc(1000, 5000, 0.01)