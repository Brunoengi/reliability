import numpy as np
from main import *

## Expected β = 2.74
def gfunction(x, d):
    angle_deg = 45.0 - x[1] / 2.0
    tan2 = np.tan(np.radians(angle_deg)) ** 2   # tan^2(45 - X2/2),  converting to radians
    g = 27.668 * x[0] + 18.595 * x[2] - 121.5 * x[0] * tan2
    return g

xvar = [
    {'varname': 'x1', 'vardist': 'normal', 'varmean': 16.0, 'varstd': 1.12},
    {'varname': 'x2', 'vardist': 'normal', 'varmean': 30.0, 'varstd': 3.0},
    {'varname': 'x3', 'vardist': 'normal', 'varmean': 25.0, 'varstd': 1.0},
]

dvar = [{'varname': 'factor1', 'varvalue': 1.00}]

example5 = Reliability(xvar, dvar, gfunction)
example5.mc(1000, 5000, 0.01)