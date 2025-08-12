from main import *

## Expected β ≈ 2.13 (≈ -> less reliable answer)

def gfunction(x, d):
    sum_sq = sum([xi**2 for xi in x])
    g = 2 + 0.015 * sum_sq - x[9]  # x[9] is X_10
    return g

xvar = [
    {'varname': f'x{i+1}', 'vardist': 'normal', 'varmean': 0.0, 'varstd': 1.0}
    for i in range(10)
]

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

example = Reliability(xvar, dvar, gfunction)
example.mc(1000, 5000, 0.01)