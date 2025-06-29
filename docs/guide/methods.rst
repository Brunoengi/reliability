Reliability Methods
================================

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
    reliability_problem = Reliability(xvar, dvar, gfunction)

It is now possible to choose a method contained in the Reliability class to solve the reliability problem. 


Programmed Reliability Methods
-------------------------------------------------

FORM (First-order reliability method): 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Algorithm FORM-iHLRF: Normal equivalente transformation**

.. code-block:: bash

  .form(iHLRF=True, tolerance)

**Algorithm FORM-iHLRF: Direct mapping to standard Gaussian space**

.. code-block:: bash

  .form2(iHLRF=True, tolerance)


References: 

HASOFER, A. M.; LIND, N. C. Exact and invariant second-moment code format. Journal of the Engineering Mechanics Division, v. 100, n. 1, p. 111–121, 1974.

DITLEVSEN, O.; MADSEN, H. O. Structural reliability methods. Chichester: Wiley, 1996.

.. raw:: html

   <br>

MCS (Monte Carlo Simulation): 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Monte Carlo Brute Force = no adaptive technique**

.. code-block:: bash

  .mc(nc, ns, delta_lim)

References: 

METROPOLIS, N.; ULAM, S. The Monte Carlo Method. Journal of the American Statistical Association, v. 44, n. 247, p. 335–341, 1949.

RUBINSTEIN, R. Y.; KROESE, D. P. Simulation and the Monte Carlo Method. 3. ed. Hoboken: Wiley, 2016.

.. raw:: html

   <br>

MCS (Monte Carlo Simulation - Variance Reduction Techniques): 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Importance sampling based on project point**

.. code-block:: bash

  .sampling_project_point(nc, ns, delta_lim)

Reference: 

BORGUND, U.; BUCHER, C. G. Importance sampling procedure using design point – ISPUD: user’s manual. Innsbruck: Institut für Mechanik, Universität Innsbruck, 1986.

.. raw:: html

   <br>

**Importance sampling with adaptive technique - Search-based importance sampling**

.. code-block:: bash

  .adaptive(nc, ns, delta_lim)

Reference: Melchers, R.E. Search-based importance sampling. Structural Safety, 9 (1990) 117-128

.. raw:: html

   <br>

**Importance sampling with adaptive technique - Iterative procedure**

.. code-block:: bash

  .bucher(nc, ns, delta_lim)

Reference: BUCHER, C.G. Adaptive sampling – an iterative fast Monte Carlo procedure. Structural safety, v. 5, n. 2, p. 119-126, 1988.

.. raw:: html

   <br>


**Enhanced Sampling** 

.. code-block:: bash

  .sampling_enhanced(nc, ns, delta_lim)

Reference: Naess A, Leira BJ, Batsevych O, 2009: System reliability analysis by enhanced Monte Carlo simulation, Structural Safety 31, 349-355.


List of parameters:
--------------------------------------------

.. raw:: html

   <br>


=========================  =========================  =========================
Parameter                        Type                      Recomendation
=========================  =========================  =========================
nc                            integer                       50 ≥ nc ≥ 200                  
ns                            integer                     2000 ≥ ns ≥ 10000 
delta_lim                     float                       0.005 ≥ delta_lim ≥0.05
tolerance                     float                       1e-6
=========================  =========================  =========================
