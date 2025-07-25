import scipy.optimize
import numpy as np
import scipy.linalg
from scipy.stats import norm, multivariate_normal



class RandomVariablesGenerator:
  def __init__(self, parent):
      ## Get all properties about Monte Carlo Methods
      self.parent = parent

      ## Get all properties about Reliability
      self.reliability = parent.reliability

  def main(self, ns):
        """
        Method to generate random variables, check correlation matrix to use var_gen (correlated variables) or var_rvs (uncorrelated variables )

        """
        # Number of variables of the problem
        ns = int(ns)
        
        # Get index of correlated and uncorrelated variables
        index_correlated, index_uncorrelated = self.reliability.correlation.correlation_summary()
        total_index = len(index_correlated) + len(index_uncorrelated) 

        # Empty matrix 
        x = np.empty((ns, total_index))
        wpc = np.ones(ns)
        fxc = np.ones(ns)
        wpu = np.ones(ns)
        fxu = np.ones(ns)
        
        if index_correlated:
          xpc, wpc, fxc = self.var_gen(ns, index_correlated)
        
        if index_uncorrelated:
          xpu, wpu, fxu = self.var_rvs(ns, index_uncorrelated)
        
        for i, idx in enumerate(index_correlated):
          x[:, idx] = xpc[:, i]

        # Fill in the uncorrelated columns
        if index_correlated:
          for i, idx in enumerate(index_correlated):
            x[:, idx] = xpc[:, i]

        if index_uncorrelated:
          for i, idx in enumerate(index_uncorrelated):
            x[:, idx] = xpu[:, i]

        wp = wpc * wpu
        fx = fxc * fxu

        return x, wp, fx

  def var_gen(self, ns, indexes_correlated_xvar, iprint=False):
      """
      Random variables generator for Monte Carlo Simulation methods, only to correlated variables
      """
      xvar_correlated = [self.reliability.xvar[i] for i in indexes_correlated_xvar]
      nxvar_correlated = len(xvar_correlated)

      # Generate initial matrices and variables
      uk_cycle = np.random.rand(ns, nxvar_correlated)
      x = np.zeros((ns, nxvar_correlated))
      weight = np.ones(ns)
      fxixj = np.ones(ns)
      zf = np.zeros((ns, nxvar_correlated))

      # Correlation submatrix for correlated variables
      matrix = self.reliability.correlation.Rz_rectify[np.ix_(indexes_correlated_xvar, indexes_correlated_xvar)]

      # Cholesky to generate correlated Gaussian samples
      L = scipy.linalg.cholesky(matrix, lower=True)
      yk = norm.ppf(uk_cycle)       # Independent standard normals
      zk = yk @ L.T                 # Correlated normals

      for i, distribution in enumerate(xvar_correlated):
          zk_col = zk[:, i]
          x[:, i], fx, hx, zf[:, i] = distribution.transform(zk_col)

          # Update weights and fxixj for variable i
          w, fx_over_norm = distribution.update_weights(fx, hx, zf[:, i], zk_col)
          weight *= w
          fxixj *= fx_over_norm

      # Final correction with multivariate PDFs
      norm_multivar = multivariate_normal(cov=matrix)
      phif = norm_multivar.pdf(zf)
      phih = norm_multivar.pdf(zk)

      weight *= phif / phih
      fxixj *= phif

      return x, weight, fxixj 
  
  def var_rvs(self, ns, indexes_uncorrelated_xvar, iprint=False):
    """
    Random variables generator for Monte Carlo Simulation methods (uncorrelated variables only)
    """

    xvar_uncorrelated = [self.reliability.xvar[i] for i in indexes_uncorrelated_xvar]
    nxvar_uncorrelated = len(xvar_uncorrelated)

    x = np.zeros((ns, nxvar_uncorrelated))
    weight = np.ones(ns)
    fxixj = np.ones(ns)

    for i, distribution in enumerate(xvar_uncorrelated):
            x[:, i], fx, hx = distribution.sample_direct(ns)
            w = fx / hx
            weight *= w
            fxixj *= fx
        
    return x, weight, fxixj
  

