from main import *


def gfunction(x, d):

    g = d[0] * x[0] - (d[1] * x[1])
    return g

# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'R', 'vardist': 'uniform', 'parameter1': 0, 'parameter2': 10},
    {'varname': 'S', 'vardist': 'uniform', 'parameter1': 0, 'parameter2': 5},
    
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
