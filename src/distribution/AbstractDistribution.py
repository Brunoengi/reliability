from abc import ABC, abstractmethod

from utils.validate.domain_types.validate_xvar import ValidateXvar

class AbstractDistribution(ABC):
  def __init__ (self, props: dict):
      
    ValidateXvar(props)
  
    self.set_properties(props)
    self.set_initial_values()
    

  @abstractmethod
  def validate_specific_parameters(self, props):
    pass

  def set_properties(self, props):
    for key, value in props.items():
      setattr(self, key, value)
      
  # @abstractmethod   
  # def instrumental_properties(self, varhmean):
  #   pass
  
  def set_initial_values(self, nsigma = 1):
      # Checks whether 'varhmean' was provided; if not, use 'varmean'
      self.varhmean = float(getattr(self, 'varhmean', self.varmean))
  
      # Tests whether 'varstd' exists as an attribute and calculates 'varcov' if possible
      if hasattr(self, 'varstd'):
        self.varcov = float(self.varstd / self.varmean) if self.varmean > 0 else 1.00
      else:
        self.varstd = float(self.varcov * self.varmean)

      self.nsigma = nsigma
      self.mufx = self.varmean
      self.muhx = self.varhmean
      self.sigmafx = self.varstd
      self.sigmahx = nsigma * self.sigmafx
      
      