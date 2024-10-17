# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 12:32:15 2024

@author: BrunoTeixeira
"""

from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.Dictionary import ValidateDictionary
from utils.validate.Class import ValidateClass

class Normal(AbstractDistribution):
  def __init__(self, dictionaryInfo):
      
    ValidateDictionary.is_dictionary(dictionaryInfo)
    ValidateDictionary.has_keys(dictionaryInfo,'varmean')
    ValidateDictionary.is_float(dictionaryInfo, 'varmean')
    ValidateDictionary.check_keys_count(dictionaryInfo, 1, 'varcov', 'varstd')
    super().__init__(dictionaryInfo)
    ValidateClass.has_invalid_key(self,'varname','vardist','varmean','varcov','varstd','varhmean')
      

# teste1 = Normal({'varname': 'fc', 'vardist': 'normal', 'varmean': 23.00, 'varcov': 0.15})
# teste2 = Normal({'varname': 'b', 'vardist': 'normal', 'varmean': 0.20, 'varstd': 0.012})
# teste3 = Normal({'varname': 'X1', 'vardist': 'normal', 'varmean': 0.00, 'varstd': 1.00, 'varhmean': +1.9149})

# for chave, valor in teste1.__dict__.items():
#     print(f"{chave}: {valor}")

# for chave, valor in teste2.__dict__.items():
#     print(f"{chave}: {valor}")

# for chave, valor in teste3.__dict__.items():
#     print(f"{chave}: {valor}")