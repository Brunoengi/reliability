# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 22:11:47 2024

@author: BrunoTeixeira
"""
from utils.validate.Dictionary import ValidateDictionary

class AbstractDistribution():
    def __init__ (self, dictionaryInfo):
        
      ValidateDictionary.has_keys(dictionaryInfo, 'varname','vardist')
  
      self.set_properties(dictionaryInfo)
      self.set_initial_values()

    def set_properties(self, dictionaryInfo):
      for key, value in dictionaryInfo.items():
        setattr(self, key, value)
    
    def set_initial_values(self):
        # Checks whether 'varhmean' was provided; if not, use 'varmean'
        self.varhmean = float(getattr(self, 'varhmean', self.varmean))
    
        # Tests whether 'varstd' exists as an attribute and calculates 'varcov' if possible
        if hasattr(self, 'varstd'):
            self.varcov = float(self.varstd / self.varmean) if self.varmean > 0 else 1.00
        else:
            self.varstd = float(self.varcov * self.varmean)
