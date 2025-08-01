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
      
    self.validate_specific_parameters(props)
    super().__init__(props)
    ValidateClass.has_invalid_key(self,'varname','vardist','varmean','varcov','varstd','varhmean')

  def validate_specific_parameters(self, props):
    ValidateDictionary.is_dictionary(props)
    ValidateDictionary.has_keys(props,'varmean')
    ValidateDictionary.is_float(props, 'varmean')
    ValidateDictionary.check_keys_count(props, 1, 'varcov', 'varstd')
    ValidateDictionary.check_if_exists(props, 'varcov', lambda d, k: ValidateDictionary.is_greater_or_equal_than(d, k, 0))
      

