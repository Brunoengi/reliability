from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.base_types.validate_dictionary import ValidateDictionary
from utils.validate.base_types.validate_class import ValidateClass 
from math import sqrt, pi
from scipy.stats import norm, gumbel_r
import numpy as np

class Gumbel(AbstractDistribution):
  def __init__(self, props: dict):
    
    ##self.validate_specific_parameters(props)
    super().__init__(props)
    
    self.alphafn = pi / sqrt(6) / self.sigmafx
    self.euler_gamma = 0.5772156649015329
    self.ufn = self.mufx - self.euler_gamma / self.alphafn
    self.betafn = 1 / self.alphafn
    self.alphahn = pi / sqrt(6) / self.sigmahx
    self.uhn = self.muhx - self.euler_gamma / self.alphahn
    self.betahn = 1 / self.alphahn
    
    
    ##ValidateClass.has_invalid_key(self, 'varname', 'vardist', 'varmean', 'varstd', 'varcov', 'varhmean')

  def validate_specific_parameters(self, props):
    ValidateDictionary.is_dictionary(props)
    ValidateDictionary.check_possible_arrays_keys(props, ['varmean', 'varstd'], ['varmean', 'varcov'])
    ValidateDictionary.check_if_exists(props, 'varcov', lambda d, k: ValidateDictionary.is_greater_or_equal_than(d, k, 0))
    
  def instrumental_properties(self, varhmean):
    self.varhmean = varhmean
    self.muhx = self.varhmean
    self.uhn = self.muhx - self.euler_gamma / self.alphahn
    
  def x_uncorrelated(self, ns):
    return gumbel_r.rvs(loc=self.uhn, scale=self.betahn, size=ns)
  
  def x_correlated(self, zk_col):
    uk = norm.cdf(zk_col)
    return self.uhn - self.betahn * np.log(np.log(1 / uk))
  
  def fx_uncorrelated(self, x):
    return gumbel_r.pdf(x, self.ufn, self.betafn)
  
  def fx_correlated(self, x):
    return self.fx_uncorrelated(x)
  
  def hx_uncorrelated(self, x):
    return gumbel_r.pdf(x, self.uhn, self.betahn)
  
  def hx_correlated(self, x):
    return self.hx_uncorrelated(x)
  
  def zf_correlated(self, x):
    cdfx = gumbel_r.cdf(x, self.ufn, self.betafn)
    return norm.ppf(cdfx)