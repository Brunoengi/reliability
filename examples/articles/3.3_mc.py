## Expected β = 2.658
from main import *

def gfunction(x, d):

    g = d[0]*x[0]*x[1]-d[1]*x[2]
    return g

xvar = [
    {'varname': 'Y', 'vardist': 'lognormal', 'varmean': 40.00, 'varcov': 0.125},
    {'varname': 'Z', 'vardist': 'lognormal', 'varmean': 50.00, 'varcov': 0.05},
    {'varname': 'M', 'vardist': 'gumbel', 'varmean': 1000.00, 'varcov': 0.20}
]

dvar = [
    {'varname': 'gamma1', 'varvalue': 1.00},
    {'varname': 'gamma2', 'varvalue': 1.00}
]

corrmatrix = [[1.00, 0.40, 0.00],
              [0.40, 1.00, 0.00],
              [0.00, 0.00, 1.00]]


beam = Reliability(xvar, dvar, gfunction, None, corrmatrix)
beam.mc(100, 10000, 0.005)
