from struct_reliability import *


def gfunction(x, d):

    g = d[0] * x[0] - (d[1] * x[1])
    return g

# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'R', 'vardist': 'uniform', 'varmean': 5.2, 'varstd': 0.1},
    {'varname': 'S', 'vardist': 'uniform', 'varmean': 5, 'varstd': 0.1},
    
]
# Design variables

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
    {'varname': 'factor2', 'varvalue': 1.00},

]
#
# MCS method
#
storm = Reliability(xvar, dvar, gfunction)
storm.mc(100, 5000, 0.005)
