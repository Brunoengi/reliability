from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.base_types.validate_dictionary import ValidateDictionary
from utils.validate.base_types.validate_class import ValidateClass 
import numpy as np
from scipy.stats import uniform, norm
from scipy.optimize import fsolve

class Uniform(AbstractDistribution):
  def __init__(self, props):

    if(set(['parameter1','parameter2']).issubset(props)):
      self.a = props['parameter1']
      self.b = props['parameter2']
    
      self.varmean = float((self.a + self.b) / 2)
      self.varstd = float((self.b - self.a) / np.sqrt(12))
      
      ##Tem que ver isso se não está sobreescrevendo
      self.varhmean = float(self.varmean)
      
    elif(set(['varmean','varstd']).issubset(props)):
      self.varmean = props['varmean']
      self.varstd = props['varstd']
      
      delta = np.sqrt(3) * self.varstd
      self.a = self.varmean - delta
      self.b = self.varmean + delta
    
    else:
      raise ValueError("Uniform distribution requires either (parameter1, parameter2) or (varmean, varstd)")
    
    super().__init__(props)
    self.ah, self.bh = fsolve(self.uniform_limits, (1, 1), args= (self.muhx, self.sigmahx)) 
    ##ValidateClass.has_invalid_key(self, 'varname', 'vardist', 'varmean', 'varcov', 'varstd', 'varhmean', 'parameter1', 'parameter2')

  def validate_specific_parameters(self, props):
    ValidateDictionary.check_possible_arrays_keys(props, ['varmean','varstd'], ['parameter1', 'parameter2'])
    
  def instrumental_properties(self, varhmean):
    self.varhmean = varhmean
    self.muhx = self.varhmean
    self.ah, self.bh = fsolve(self.uniform_limits, (1, 1), args= (self.muhx, self.sigmahx))
    
  def x_uncorrelated(self, ns):
    return uniform.rvs(loc=self.ah, scale= (self.bh-self.ah), size = ns)
  
  def x_correlated(self, zk_col):
    self.uk = norm.cdf(zk_col)
    return self.ah + (self.bh - self.ah) * self.uk
  
  def fx_uncorrelated(self, x):
    return uniform.pdf(x, self.a, self.b-self.a)
  
  def fx_correlated(self, x):
    return self.fx_uncorrelated(x)
  
  def hx_uncorrelated(self, x):
    return uniform.pdf(x, self.ah, self.bh-self.ah)
  
  def hx_correlated(self, x):
    return self.hx_uncorrelated(x)
  
  def zf_correlated(self, x):
    return norm.ppf(self.uk)

  @staticmethod
  def uniform_limits(vars, mux, sigmax):
    a, b = vars
    eq1 = (a + b) / 2 - mux
    eq2 = (b - a) / np.sqrt(12.) - sigmax
    return [eq1, eq2]