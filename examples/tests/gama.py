from struct_reliability import *
import math

def gfunction(x, d):

    g = d[0] * x[0] - (d[1] * x[1])
    return g

# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'R', 'vardist': 'gamma', 'varmean': 11/22.22, 'varstd': math.sqrt(11/22.22**2)},
    {'varname': 'S', 'vardist': 'gamma', 'varmean': 11/22.22, 'varstd': math.sqrt(11/22.22**2)},
    
]
# Design variables

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
    {'varname': 'factor2', 'varvalue': 1.00},

]

corrmatrix = [[1.00, 0.75],
              [0.75, 1.00]]

# MCS method
#
column = Reliability(xvar, dvar, gfunction, None, corrmatrix)
column.mc(100, 5000, 0.005)
