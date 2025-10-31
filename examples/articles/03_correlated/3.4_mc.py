## Expected β = 3.397
from struct_reliability import *

def gfunction(x, d):

    g = d[0]*x[0]*x[1]-d[1]*x[2]
    g= 0.02 - 8 * (x[0]/(x[1] * x[2]))
    return g

xvar = [
    {'varname': 'q', 'vardist': 'lognormal', 'varmean': 1000, 'varcov':0.2},
    {'varname': 'E', 'vardist': 'lognormal', 'varmean': 2e10, 'varcov': 0.05},
    {'varname': 'I', 'vardist': 'lognormal', 'varmean': 3.9025e-5, 'varcov': 0.1}
]

dvar = [
    {'varname': 'gamma1', 'varvalue': 1.00},
    {'varname': 'gamma2', 'varvalue': 1.00}
]

corrmatrix = [[1.00, 0.20, 0.20],
              [0.20, 1.00, 0.20],
              [0.20, 0.20, 1.00]]


beam = Reliability(xvar, dvar, gfunction, None, corrmatrix)
beam.mc(10000, 5000, 0.01)
