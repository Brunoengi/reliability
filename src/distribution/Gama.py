from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.base_types.validate_dictionary import ValidateDictionary
from utils.validate.base_types.validate_class import ValidateClass 
import scipy.optimize
import numpy as np
import scipy.linalg
from math import log
from scipy.stats import norm, uniform, lognorm, gumbel_r, invweibull, weibull_min, beta as beta_dist, gamma as gamma_dist, multivariate_normal
from scipy.optimize import fsolve, newton
from scipy.linalg import cholesky
from math import sqrt, pi, log
from scipy.special import gamma

class Gama(AbstractDistribution):
  def __init__(self, props: dict):

    ##self.validate_specific_parameters(props)
    super().__init__(props)
    self.deltafx = self.sigmafx / self.mufx
    self.deltahx = self.sigmahx / self.muhx
    self.k = 1 / self.deltafx ** 2
    self.v = self.k / self.mufx
    self.a, self.loc, self.scale = self.k, 0, 1 / self.v
    self.kh = 1 / self.deltahx ** 2
    self.vh = self.kh / self.muhx
    self.ah, self.loch, self.scaleh = self.kh, 0, 1 / self.vh

  def validate_specific_parameters(self, props):
    ValidateDictionary.is_dictionary(props)
    ValidateDictionary.check_possible_arrays_keys(props, ['varmean', 'varcov'], ['varmean', 'varstd'])
    ValidateDictionary.check_if_exists(props, 'varcov', lambda d, k: ValidateDictionary.is_greater_or_equal_than(d, k, 0))
    
  def instrumental_properties(self, varhmean):
    self.varhmean = varhmean
    self.muhx = self.varhmean
    self.deltahx = self.sigmahx / self.muhx
    self.kh = 1 / self.deltahx ** 2
    self.vh = self.kh / self.muhx
    self.ah, self.loch, self.scaleh = self.kh, 0, 1 / self.vh
  
  def x_correlated(self, zk_col):
    uk = norm.cdf(zk_col)
    return gamma_dist.ppf(uk, self.ah, loc=self.loch, scale=self.scaleh)
  
  def fx_correlated(self, x):
    return gamma_dist.pdf(x, self.a, loc=self.loc, scale=self.scale)
  
  def hx_correlated(self, x):
    return gamma_dist.pdf(x, self.ah, loc=self.loch, scale=self.scaleh)
  
  def zf_correlated(self, x):
    cdfx = gamma_dist.cdf(x, self.a, loc=self.loc, scale=self.scale)
    return norm.ppf(cdfx)
  
  def x_uncorrelated(self, ns):
    return gamma_dist.rvs(self.ah, self.loch, self.scaleh, size=ns)
  
  def fx_uncorrelated(self,x):
    return gamma_dist.pdf(x, self.a, self.loc, self.scale)
  
  def hx_uncorrelated(self, x):
    return gamma_dist.pdf(x, self.ah, self.loch, self.scaleh)
  
