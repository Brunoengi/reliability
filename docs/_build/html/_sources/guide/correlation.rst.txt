Correlation between variables
=====================================


To solve problems in which there are correlated variables, it is possible to insert an correlation matrix containing the correlation coefficients between pairs of pairs of analyzed variables.

To construct the matrix, some rules must be followed:

1 - the main diagonal must be 1, considering that the display of an identified variable must be 1;

2 - the other coefficients of the matrix must be values ​​greater than or equal to -1 and less than or equal to 1;

3 - the number of rows and columns of the matrix must be equal to the number of random variables of the reliability problem that you want to solve.


For example, you can define a situation with 4 random variables, as shown below:

.. code-block:: bash

  xvar = [
      {'varname': 'R', 'vardist': 'normal', 'varmean': 975.00, 'varcov': 0.15},
      {'varname': 'G', 'vardist': 'normal', 'varmean': 200.00, 'varcov': 0.07},
      {'varname': 'Q', 'vardist': 'normal', 'varmean': 300.00, 'varcov': 0.12},
      {'varname': 'w', 'vardist': 'normal', 'varmean': 150.00, 'varcov': 0.20}
  ]


In this case, the correlation matrix must contain 4 rows and 4 columns, since the number of calculated variables is equal to 4.
An example of brightness between variables is shown below:

.. math::

   \rho_{x_{12}} = \rho_{x_{21}} = 0{,}8

.. math::

   \rho_{x_{23}} = \rho_{x_{32}} = 0{,}3

Considering:

:math:`x_{1} = R \hspace{2em} x_{2} = G \hspace{2em} x_{3} = Q \hspace{2em} x_{4} = W`


In this context, it is possible to define the correlation matrix :math:`R_{x}`:

:math:`R_x = \begin{bmatrix}
1{,}0 & 0{,}8 & 0{,}0 & 0{,}0 \\
0{,}8 & 1{,}0 & 0{,}3 & 0{,}0 \\
0{,}0 & 0{,}3 & 1{,}0 & 0{,}0 \\
0{,}0 & 0{,}0 & 0{,}0 & 1{,}0
\end{bmatrix}`

The complete problem is presented below:


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