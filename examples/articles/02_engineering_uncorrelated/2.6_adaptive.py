from main import *

## Expected β = 4.78
def gfunction(x, d):
    numerator = np.sqrt(3 * (1 - 0.3**2))
    denom = np.pi * x[0] * (x[1]**2) * (np.cos(x[2])**2)
    factor = numerator / denom
    term = x[5] / 0.66 + x[4] / (0.41 * x[3]) 
    g = 1 - factor * term
    return g

xvar = [
    {'varname': 'x1', 'vardist': 'normal', 'varmean': 7.0e10, 'varstd': 3.50e9},
    {'varname': 'x2', 'vardist': 'normal', 'varmean': 2.5e-3, 'varstd': 1.25e-4},
    {'varname': 'x3', 'vardist': 'normal', 'varmean': 0.524, 'varstd': 0.01048},
    {'varname': 'x4', 'vardist': 'normal', 'varmean': 0.90, 'varstd': 0.0225},
    {'varname': 'x5', 'vardist': 'normal', 'varmean': 8e4, 'varstd': 6.4e3},
    {'varname': 'x6', 'vardist': 'normal', 'varmean': 7e4, 'varstd': 5.6e3},
]

dvar = [{'varname': 'factor1', 'varvalue': 1.00}]

example6 = Reliability(xvar, dvar, gfunction)
example6.adaptive(1000, 5000, 0.01)