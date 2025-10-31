from struct_reliability import *

## Expected β = 2.59
def gfunction(x, d):
    g = 0.03 - (((x[0] * (x[1]**2)) / 2) * ((3.81/(x[3]*x[5])) + (1.13/(x[2]*x[4]))))
    return g

xvar = [
    {'varname': 'x1', 'vardist': 'normal', 'varmean': 2.0e4, 'varstd': 1.4e3},
    {'varname': 'x2', 'vardist': 'normal', 'varmean': 12, 'varstd': 0.12},
    {'varname': 'x3', 'vardist': 'normal', 'varmean': 9.82e-4, 'varstd': 5.9825e-5},
    {'varname': 'x4', 'vardist': 'normal', 'varmean': 0.04, 'varstd': 4.8e-3},
    {'varname': 'x5', 'vardist': 'normal', 'varmean': 1.0e11, 'varstd': 1.0e9},
    {'varname': 'x6', 'vardist': 'normal', 'varmean': 2.0e10, 'varstd': 1.2e9},
]

dvar = [{'varname': 'factor1', 'varvalue': 1.00}]

example7 = Reliability(xvar, dvar, gfunction)
example7.mc(100000, 5000, 0.01)