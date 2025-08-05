from scipy.stats import norm

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
      
  def validate_specific_parameters(self, props):
    ValidateDictionary.is_dictionary(props)
    ValidateDictionary.has_keys(props,'varmean')
    ValidateDictionary.is_float(props, 'varmean')
    ValidateDictionary.check_if_exists(props, 'varcov', lambda d, k: ValidateDictionary.is_greater_or_equal_than(d, k, 0))
      
  def x_correlated(self, zk_col):
    return self.muhx + self.sigmahx * zk_col
  
  def x_uncorrelated(self, ns):
    return norm.rvs(loc=self.muhx, scale=self.sigmahx, size=ns)
  
  def fx(self, x):
    return norm.pdf(x, self.mufx, self.sigmafx)
  
  def hx(self, x):
    return norm.pdf(x, self.muhx, self.sigmahx)
  
  def zf(self, x):
    return (x - self.mufx) / self.sigmafx
