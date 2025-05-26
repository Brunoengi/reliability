Reliability Methods
=========================

The last step in developing a reliability problem is to instantiate a class and choose a resolution method.
To instantiate the class, 3 mandatory parameters are required and there is the possibility of two optional parameters.
The mandatory parameters are precisely the set of random variables (:doc:`variables`), the design variables (:doc:`design`) and the limit state function (:doc:`limit`).

An example containing the basic situation is described below:

  .. code-block:: bash

    from main import Reliability

    def gfunction(x, d):

      g = d[0]*x[0]*x[1]-d[1]*x[2]
      return g


    #
    # Data input
    #
    # Random variables: name, probability distribution, mean and coefficient of variation


    xvar = [
        {'varname': 'Y', 'vardist': 'normal', 'varmean': 40.00, 'varcov': 0.125},
        {'varname': 'Z', 'vardist': 'normal', 'varmean': 50.00, 'varcov': 0.05},
        {'varname': 'M', 'vardist': 'normal', 'varmean': 1000.00, 'varcov': 0.20}
    ]

    # Design variables

    dvar = [
        {'varname': 'gamma1', 'varvalue': 1.00},
        {'varname': 'gamma2', 'varvalue': 1.00}
    ]

    #
    # Instantiating the class
    #
    reliabilityProblem = Reliability(xvar, dvar, gfunction)

It is now possible to choose a method contained in the Reliability class to solve the reliability problem. 


Programmed Reliability Methods
--------------------------------------------

FORM (First-order reliability method): 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bibliography: 


Computational Methods: 


``.form``: Algorithm FORM-iHLRF. Normal equivalente transformation

``.form2``: Algorithm FORM-iHLRF. Direct mapping to standard Gaussian space

MCS (Monte Carlo Simulation): 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bibliography:

Computational Methods: 

``.mc``: Monte Carlo simulation method without adaptive technique using proprietary algorithm for variable generation

``.mc2``: Monte Carlo simulation method without adaptive technique using native Python variable generation

``.adaptive``: Monte Carlo simulations method with Importance Sampling (MC-IS) - Importance sampling with adaptive technique using proprietary algorithm for variable generation

``.adaptive2``: Monte Carlo simulations method with Importance Sampling (MC-IS) - Importance sampling with adaptive technique using native Python variable generation

``.bucher``: Monte Carlo Simulations with Importance Sampling (MC-IS) - Importance sampling with adaptive technique using proprietary algorithm for variable generation

``.bucher2``: Monte Carlo Simulations with Importance Sampling (MC-IS) - Importance sampling with adaptive technique using native Python variable generation

