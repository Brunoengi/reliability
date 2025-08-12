from main import *

## Expected β = 2.709

def gfunction(x, d):
    x1 = x[0]
    x2 = x[1]
    g = np.exp(0.4 * (x1 + 2) + 6.2) - np.exp(0.3 * x2 + 5) - 200
    return g

xvar = [
    {'varname': 'x1', 'vardist': 'normal', 'varmean': 0.0, 'varstd': 1.0},
    {'varname': 'x2', 'vardist': 'normal', 'varmean': 0.0, 'varstd': 1.0},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

example = Reliability(xvar, dvar, gfunction)
example.mc(1000, 5000, 0.01)