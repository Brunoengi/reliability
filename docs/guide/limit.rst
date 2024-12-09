Limit State Function
=========================

To define the limit state function it is necessary to create a function with two parameters, the first parameter refers to random variables, while the second parameter refers to design variables. 
In this context, it is interesting that you have already read the documentation regarding random variables and design variables.

Let's say we want to create a problem in which the limit state function is dependent on three random variables and two design variables. In this case, the first parameter, referring to the random variables, will be a list with 3 positions (from zero to two) while the parameter referring to the design variables will be a list with two positions (from zero to one).

.. code-block:: bash

    def gfunction(x, d):

      g = d[0] * x[0] * x[1] - d[1] * x[2]
      return g

In this case, x references researched variables while d references design variables, note that the variable ``x`` will be a list with three positions while the variable ``d`` will be a list with two positions.

For this situation to make sense, it is necessary to have created 3 random variables and two design variables. A possible situation is as described below:

.. code-block:: bash

  # Random variables
  xvar = [
    {'varname': 'Y', 'vardist': 'lognormal', 'varmean': 38.00, 'varcov': 0.10},
    {'varname': 'Z', 'vardist': 'normal', 'varmean': 60.00, 'varcov': 0.05},
    {'varname': 'M', 'vardist': 'frechet', 'varmean': 1000.00, 'varcov': 0.30}
  ]

  # Design variables
  dvar = [
    {'varname': 'gamma1', 'varvalue': 1.00},
    {'varname': 'gamma2', 'varvalue': 1.00}
  ]

In this case, ``d[0]`` refers to ``gamma1``, ``d[1]`` refers to ``gamma2``, ``x[0]`` refers to ``Y``, ``x[1]`` refers to ``Z``, and ``x[2]`` refers to ``M``.

Note that the limit state function needs two parameters that are lists and the return of the function will be of type float.
