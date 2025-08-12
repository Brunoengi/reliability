from main import *

## Expected β ≈ 2.29 (≈ -> less reliable answer)
def gfunction_case8(x, d):
    X1 = float(x[0])
    X2 = float(x[1])

    a = 1 - 0.01 * (1 / X1)**2 - (1 / X1)**2 - (1 / X2)**2 + (1 / (X1 * X2))**2

    invX1 = 1.0 / X1
    invX2 = 1.0 / X2
    

    numer = abs(1 - (invX2)**2)
    denom = ((a)**2 + 4 * (0.01**2) * ((invX1**2 - (1/(X1*X2**2)))**2))**0.5

    g = 27 - numer / denom
    return g

xvar_case8 = [
    {'varname': 'x1', 'vardist': 'normal', 'varmean': 1.0, 'varstd': 0.025},
    {'varname': 'x2', 'vardist': 'normal', 'varmean': 1.0, 'varstd': 0.025},
]

dvar_case8 = [{'varname': 'factor1', 'varvalue': 1.00}]

example8 = Reliability(xvar_case8, dvar_case8, gfunction_case8)
example8.mc(100000, 5000, 0.01)