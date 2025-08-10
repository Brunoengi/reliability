import scipy.optimize
import numpy as np
import scipy.linalg
from math import log
from scipy.stats import norm, uniform, lognorm, gumbel_r, invweibull, weibull_min, beta as beta_dist, gamma as gamma_dist, multivariate_normal
from scipy.optimize import fsolve, newton
from scipy.linalg import cholesky
from math import sqrt, pi, log
from scipy.special import gamma


class RandomVariablesGenerator:
  def __init__(self, parent):
      ## Get all properties about Monte Carlo Methods
      self.parent = parent

      ## Get all properties about Reliability
      self.reliability = parent.reliability

  def main(self, ns, nsigma=1.00, iprint=False):
        """
        Method to generate random variables, check correlation matrix to use var_gen (correlated variables) or var_rvs (uncorrelated variables)

        """
        # Number of variables of the problem
        ns = int(ns)
        
        # Get index of correlated and uncorrelated variables
        index_correlated, index_uncorrelated = self.reliability.correlation.correlation_summary()
        total_index = len(index_correlated) + len(index_uncorrelated) 

        # Empty matrix 
        x = np.empty((ns, total_index))

        #
        # Standard deviation multiplier for MC-IS
        
        #
        # Step 1 - Generation of the random numbers according to their appropriate distribution
        #

        wpc = np.ones(ns)
        fxc = np.ones(ns)
        wpu = np.ones(ns)
        fxu = np.ones(ns)
        
        if index_correlated:
          xpc, wpc, fxc = self.var_gen(ns, index_correlated, nsigma)
        
        if index_uncorrelated:
          xpu, wpu, fxu = self.var_rvs(ns, index_uncorrelated, nsigma)
        


        for i, idx in enumerate(index_correlated):
          x[:, idx] = xpc[:, i]

        # Preencher as colunas dos não correlacionados
        if index_correlated:
          for i, idx in enumerate(index_correlated):
            x[:, idx] = xpc[:, i]

        if index_uncorrelated:
          for i, idx in enumerate(index_uncorrelated):
            x[:, idx] = xpu[:, i]

        wp = wpc * wpu
        fx = fxc * fxu

        return x, wp, fx

  def var_gen(self, ns, indexes_correlated_xvar, nsigma=1.00, iprint=False):
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

      def update_weights(fx, hx, zf_col, zk_col):
          # Update weights and fxixj product for each variable
          standard_norm_pdf_zf = norm.pdf(zf_col)
          standard_norm_pdf_zk = norm.pdf(zk_col)
          w = (fx / standard_norm_pdf_zf) / (hx / standard_norm_pdf_zk)
                
          return w, fx / standard_norm_pdf_zf

      for i, var in enumerate(xvar_correlated):
          # Adjust std if needed
          if var['varstd'] == 0.0:
              var['varstd'] = float(var['varcov']) * float(var['varmean'])

          zk_col = zk[:, i]
          x[:, i] = self.reliability.xvarClassCorrelated[i].x_correlated(zk_col)
          fx = self.reliability.xvarClassCorrelated[i].fx_correlated(x[:, i])
          hx = self.reliability.xvarClassCorrelated[i].hx_correlated(x[:, i]) 
          zf[:, i] = self.reliability.xvarClassCorrelated[i].zf_correlated(x[:, i])

          # Update weights and fxixj for variable i
          w, fx_over_norm = update_weights(fx, hx, zf[:, i], zk_col)
          weight *= w
          fxixj *= fx_over_norm

      # Final correction with multivariate PDFs
      norm_multivarf = multivariate_normal(cov=matrix)
      phif = norm_multivarf.pdf(zf)
      norm_multivarh = multivariate_normal(cov=matrix)
      phih = norm_multivarh.pdf(zk)
      weight *= phif / phih
      fxixj *= phif

      return x, weight, fxixj 
  
  def var_rvs(self, ns, indexes_uncorrelated_xvar, nsigma=1.00, iprint=False):
    """
        Random variables generator for the Monte Carlo Simulation methods, only to uncorrelated variables
    """

    #Get only correlated variables
    xvar_uncorrelated = [self.reliability.xvar[i] for i in indexes_uncorrelated_xvar]
    nxvar_uncorrelated = len(xvar_uncorrelated)

    x = np.zeros((ns, nxvar_uncorrelated))
    weight = np.ones(ns)
    fx = np.zeros(ns)
    hx = np.zeros(ns)
    fxixj = np.ones(ns)
    
    #
    i = -1
    for var in xvar_uncorrelated:
        i += 1
        x[:, i] = self.reliability.xvarClassUncorrelated[i].x_uncorrelated(ns)
        fx = self.reliability.xvarClassUncorrelated[i].fx_uncorrelated(x[:, i])
        hx = self.reliability.xvarClassUncorrelated[i].hx_uncorrelated(x[:, i])
        weight = weight * (fx / hx)
        fxixj = fxixj * fx 
            

    return x, weight, fxixj
  

