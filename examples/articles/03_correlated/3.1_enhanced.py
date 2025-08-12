# Expected pf = 0.03522 and β = 1.809
from main import *

def gfunction(x, d):

    g = d[0] * x[0] - d[1] * x[1] 
    return g

# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'R', 'vardist': 'normal', 'varmean': 150, 'varcov': 0.1333333333333333},
    {'varname': 'S', 'vardist': 'normal', 'varmean': 120, 'varcov': 0.2083333333333333},
]
# Design variables

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
    {'varname': 'factor2', 'varvalue': 1.00},
]

# Correlation matrix

corrmatrix = [[1.00, 0.75],
              [0.75, 1.00]]

# MCS method
#
column = Reliability(xvar, dvar, gfunction, None, corrmatrix)
column.sampling_enhanced(1000, 5000, 0.01)
