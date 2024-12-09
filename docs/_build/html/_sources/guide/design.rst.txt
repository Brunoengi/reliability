Design Variables
======================

Design variables reference variables that will be treated as deterministic when analyzing the limit state function. Commonly, design variables represent reduction factors and increase coefficients to achieve minimum reliability values ​
​established in technical standards.

Design variables must be created in the form of a list of dictionaries, where each dictionary represents a design variable. There are two keys required to construct a design variable. The keys are: ``varname`` and ``varvalue``.

=========================  =========================
Key                         Type (Value)
=========================  =========================
varname                       string
varvalue                      number | float
=========================  =========================

An example of how design variables should be created is shown below. In this case, a list with 3 dictionaries will be presented, which indicates that the limit state function is dependent on 3 design variables.

  .. code-block:: bash

    dvar = [
      {'varname': 'factor1', 'varvalue': 1.00},
      {'varname': 'factor2', 'varvalue': 1.00},
      {'varname': 'factor3', 'varvalue': 1.00},
    ] 