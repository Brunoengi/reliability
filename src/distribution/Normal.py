from scipy.stats import norm, uniform

# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 12:32:15 2024

@author: BrunoTeixeira
"""

from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.base_types.validate_class import ValidateClass
from utils.validate.base_types.validate_dictionary import ValidateDictionary 

class Normal(AbstractDistribution):
  def __init__(self, props: dict):
      
    ##self.validate_specific_parameters(props)
    super().__init__(props)
    
    ##ValidateClass.has_invalid_key(self,'varname','vardist','varmean','varcov','varstd','varhmean')

  def validate_specific_parameters(self, props):
    ValidateDictionary.is_dictionary(props)
    ValidateDictionary.has_keys(props,'varmean')
    ValidateDictionary.is_float(props, 'varmean')
    ValidateDictionary.check_keys_count(props, 1, 'varcov', 'varstd')
    ValidateDictionary.check_if_exists(props, 'varcov', lambda d, k: ValidateDictionary.is_greater_or_equal_than(d, k, 0))
      
  def instrumental_properties(self, varhmean):
    self.varhmean = varhmean
    self.muhx = self.varhmean
  
  def x_uncorrelated(self, ns):
    return norm.rvs(loc=self.muhx, scale=self.sigmahx, size=ns)
  
  def x_correlated(self, zk_col):
    return self.muhx + self.sigmahx * zk_col
  
  def fx_correlated(self, x):
    return self.fx_uncorrelated(x)
  
  def fx_uncorrelated(self, x):
    return norm.pdf(x, loc=self.mufx, scale=self.sigmafx)
  
  def hx_correlated(self, x):
    return self.hx_uncorrelated(x)
  
  def hx_uncorrelated(self, x):
    return norm.pdf(x, loc=self.muhx, scale=self.sigmahx)
  
  def zf_correlated(self, x):
    return (x - self.mufx) / self.sigmafx
