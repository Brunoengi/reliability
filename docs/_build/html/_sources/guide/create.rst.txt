Simple Examples
=========================

Example 1 - Reliability in columns
--------------------------------------------

Given a reinforced concrete column subjected to the loads described below, assuming a centered compression situation.

Permanent Load (G): :math:`\mu_G = 200\,\text{kN},\hspace{1em} \sigma_G = 14\,\text{kN},\hspace{1em} \delta_G = 7\%`

Accidental Load (Q): :math:`\mu_Q = 300\,\text{kN},\hspace{1em} \sigma_Q = 36\,\text{kN},\hspace{1em} \delta_Q = 12\%`

Wind Load (W): :math:`\mu_Q = 150\,\text{kN},\hspace{1em} \sigma_Q = 30\,\text{kN},\hspace{1em} \delta_Q = 20\%`

Total Load (S): :math:`S = G + Q + W`

Resistance (R): :math:`\mu_R = 975\,\text{kN},\hspace{1em} \sigma_R = 146.25\,\text{kN},\hspace{1em} \delta_Q = 15\%`


Computational development
*********************************************

In this case, we will only use the FORM method as an example

.. code-block:: bash

  from main import Reliability

  #
  # Step 0 - Column: g(R, G, Q, W) = R-G-Q-W = 0
  #


  def gfunction(x, d):

      g = d[0] * x[0] - d[1] * x[1] - d[2] * x[2] - d[3] * x[3]
      return g



  # Data input
  #
  # Random variables: name, probability distribution, mean and coefficient of variation


  xvar = [
      {'varname': 'R', 'vardist': 'normal', 'varmean': 975.00, 'varcov': 0.15},
      {'varname': 'G', 'vardist': 'normal', 'varmean': 200.00, 'varcov': 0.07},
      {'varname': 'Q', 'vardist': 'normal', 'varmean': 300.00, 'varcov': 0.12},
      {'varname': 'w', 'vardist': 'normal', 'varmean': 150.00, 'varcov': 0.20}
  ]
  # Design variables

  dvar = [
      {'varname': 'factor1', 'varvalue': 1.00},
      {'varname': 'factor2', 'varvalue': 1.00},
      {'varname': 'factor3', 'varvalue': 1.00},
      {'varname': 'factor4', 'varvalue': 1.00}
  ]

  #
  # FORM method
  #
  column = Reliability(xvar, dvar, gfunction, None, None)
  column.form(iHLRF=True, toler=1.e-6)


Example 2 - Load capacity in beams: uncorrelated variables with distribution different from normal
------------------------------------------------------------------------------------------------------------------------------------

.. image:: ../_static/images/examples/example02.png
   :alt: Descrição da imagem
   :width: 50%
   :align: center

The plastic moment (ultimate resistance capacity in the plastic regime) of a section of a steel beam can be given by: 
:math:`M_p = YZ`

Where:

Y: is the yield stress of the steel.

Z: is the plastic modulus of the cross section.

If M is the requesting moment, the performance function will be defined as:

g(X)= YZ − M

Design Parameters:

Y: lognormal distribution - :math:`\mu_Y = 40\,\text{kN/cm²},\hspace{1em} \delta_Y = 0.125\,\hspace{1em} \sigma_Y = 5\,\text{kN/cm²}`

Z: lognormal distribution - :math:`\mu_Z = 50\,\text{cm³},\hspace{1em} \delta_Z = 0.05\,\hspace{1em} \sigma_Z = 2.5\,\text{m³}`

M: Gumbel distribution - :math:`\mu_M = 1000\,\text{kN.cm},\hspace{1em} \delta_M = 0.20\,\hspace{1em} \sigma_M = 200\,\text{kN.cm}`


Computational development
*********************************************

.. code-block:: bash

  from main import Reliability

  #
  # Step 0 - Beam: g(Y, Z, M) = Y*Z-M = 0
  #


  def gfunction(x, d):

      g = d[0]*x[0]*x[1]-d[1]*x[2]
      return g


  #
  # Data input
  #
  # Random variables: name, probability distribution, mean and coefficient of variation


  xvar = [
      {'varname': 'Y', 'vardist': 'lognormal', 'varmean': 40.00, 'varcov': 0.125},
      {'varname': 'Z', 'vardist': 'lognormal', 'varmean': 50.00, 'varcov': 0.05},
      {'varname': 'M', 'vardist': 'gumbel', 'varmean': 1000.00, 'varcov': 0.20}
  ]

  # Design variables

  dvar = [
      {'varname': 'gamma1', 'varvalue': 1.00},
      {'varname': 'gamma2', 'varvalue': 1.00}
  ]
  #
  # FORM method
  #
  beam = Reliability(xvar, dvar, gfunction)
  beam.form(iHLRF=True, toler=1.e-3)
  #

Example 3 - Load capacity in beams: correlated variables with distribution different from normal
------------------------------------------------------------------------------------------------------------------------------------

The proposed problem is similar to problem 1, the difference is that the variables are correlated. The demonstration coefficients between pairs of estimates are presented below:

:math:`\rho_{x_{12}} = \rho_{x_{21}} = 0{,}8\hspace{2em} \rho_{x_{23}} = \rho_{x_{32}} = 0{,}3`

Considering:

:math:`x_{1} = R \hspace{2em} x_{2} = G \hspace{2em} x_{3} = Q \hspace{2em} x_{4} = W`


In this context, it is possible to define the correlation matrix :math:`R_{x}`:

:math:`R_x = \begin{bmatrix}
1{,}0 & 0{,}8 & 0{,}0 & 0{,}0 \\
0{,}8 & 1{,}0 & 0{,}3 & 0{,}0 \\
0{,}0 & 0{,}3 & 1{,}0 & 0{,}0 \\
0{,}0 & 0{,}0 & 0{,}0 & 1{,}0
\end{bmatrix}`

Computational development
*********************************************

.. code-block:: bash
  
  from main import Reliability

  #
  # Step 0 - Column: g(R, G, Q, W) = R-G-Q-W = 0
  #


  def gfunction(x, d):

      g = d[0] * x[0] - d[1] * x[1] - d[2] * x[2] - d[3] * x[3]
      return g


  # Data input
  #
  # Random variables: name, probability distribution, mean and coefficient of variation

  xvar = [
      {'varname': 'R', 'vardist': 'normal', 'varmean': 975.00, 'varcov': 0.15},
      {'varname': 'G', 'vardist': 'normal', 'varmean': 200.00, 'varcov': 0.07},
      {'varname': 'Q', 'vardist': 'normal', 'varmean': 300.00, 'varcov': 0.12},
      {'varname': 'w', 'vardist': 'normal', 'varmean': 150.00, 'varcov': 0.20}
  ]
  # Design variables

  dvar = [
      {'varname': 'factor1', 'varvalue': 1.00},
      {'varname': 'factor2', 'varvalue': 1.00},
      {'varname': 'factor3', 'varvalue': 1.00},
      {'varname': 'factor4', 'varvalue': 1.00}
  ]

  # Correlation matrix

  corrmatrix = [[1.00, 0.80, 0.00, 0.00],
                [0.80, 1.00, 0.30, 0.00],
                [0.00, 0.30, 1.00, 0.00],
                [0.00, 0.00, 0.00, 1.00]]

  # FORM adaptative method
  #
  column = Reliability(xvar, dvar, gfunction, None, corrmatrix)
  column.form(iHLRF=True, toler=1e-6)
  #