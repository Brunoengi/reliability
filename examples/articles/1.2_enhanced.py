import numpy as np
from main import *

## Expected β = 2.247

def gfunction(x, d):
    linear = x[0] + 2*x[1] + 2*x[2] + x[3] - 5*x[4] - 5*x[5]
    sinusoidal = 0.001 * np.sum(np.sin(100 * np.array(x[:6])))
    
    return linear + sinusoidal

xvar = [
    {'varname': 'x1', 'vardist': 'lognormal', 'varmean': 120, 'varstd': 12},
    {'varname': 'x2', 'vardist': 'lognormal', 'varmean': 120, 'varstd': 12},
    {'varname': 'x3', 'vardist': 'lognormal', 'varmean': 120, 'varstd': 12},
    {'varname': 'x4', 'vardist': 'lognormal', 'varmean': 120, 'varstd': 12},
    {'varname': 'x5', 'vardist': 'lognormal', 'varmean': 50, 'varstd': 15},
    {'varname': 'x6', 'vardist': 'lognormal', 'varmean': 40, 'varstd': 12},
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

example = Reliability(xvar, dvar, gfunction)
example.sampling_enhanced(10000, 5000, 0.01)