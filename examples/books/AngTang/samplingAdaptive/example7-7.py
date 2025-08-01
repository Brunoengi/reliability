"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.7 - Linear limit state function with non-normal independent variables
@author: MVREAL
"""
from main import *

#
# Step 0 - Column: g(R, G, Q, W) = R-G-Q-W = 0
#


def gfunction(x, d):

    g = d[0] * x[0] - d[1] * x[1] - d[2] * x[2] - d[3] * x[3]
    return g


#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'R', 'vardist': 'lognormal', 'varmean': 975.00, 'varcov': 0.15},
    {'varname': 'G', 'vardist': 'normal', 'varmean': 200.00, 'varcov': 0.07},
    {'varname': 'Q', 'vardist': 'gumbel', 'varmean': 300.00, 'varcov': 0.12},
    {'varname': 'w', 'vardist': 'gumbel', 'varmean': 150.00, 'varcov': 0.20}
]

# Design variables

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
    {'varname': 'factor2', 'varvalue': 1.00},
    {'varname': 'factor3', 'varvalue': 1.00},
    {'varname': 'factor4', 'varvalue': 1.00}
]
#
# MC-IS adaptive method
#
column = Reliability(xvar, dvar, gfunction)
column.adaptive(100, 5000, 0.03)
#