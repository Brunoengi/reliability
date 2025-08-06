from main import *


def gfunction(x, d):

    g = d[0] * x[0] - (d[1] * x[1])
    return g

# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'R', 'vardist': 'weibull', 'varmean': 10, 'varstd': 2, 'parameter3': 5},
    {'varname': 'S', 'vardist': 'weibull', 'varmean': 10, 'varstd': 2, 'parameter3': 5},
    
]
# Design variables

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
    {'varname': 'factor2', 'varvalue': 1.00},

]

# MCS method
#
column = Reliability(xvar, dvar, gfunction, None, None)
column.mc(100, 5000, 0.005)
