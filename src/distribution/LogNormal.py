from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.base_types.validate_dictionary import ValidateDictionary
from utils.validate.base_types.validate_class import ValidateClass 
from scipy.stats import lognorm
import numpy as np

class LogNormal(AbstractDistribution):
  def __init__(self, props: dict):

    ##self.validate_specific_parameters(props)
    super().__init__(props)
    self.zetafx = np.sqrt(np.log(1.00 + (self.sigmafx / self.mufx) ** 2))
    self.lambdafx = np.log(self.mufx) - 0.5 * self.zetafx ** 2
    self.zetahx = np.sqrt(np.log(1.00 + (self.sigmahx / self.muhx) ** 2))
    self.lambdahx = np.log(self.muhx) - 0.5 * self.zetahx ** 2
    ##ValidateClass.has_invalid_key(self, 'varname', 'vardist', 'varmean', 'varstd', 'varcov', 'varhmean')

  def validate_specific_parameters(self, props):
    ValidateDictionary.is_dictionary(props)
    ValidateDictionary.check_possible_arrays_keys(props,['varmean','varstd'], ['varmean','varcov'])
    ValidateDictionary.check_if_exists(props, 'varcov', lambda d, k: ValidateDictionary.is_greater_or_equal_than(d, k, 0))
    
  def x_uncorrelated(self, ns):
    return lognorm.rvs(s=self.zetahx, loc=0.00, scale=np.exp(self.lambdahx), size=ns)
  
  def x_correlated(self, zk_col):
    return np.exp(self.lambdahx + zk_col * self.zetahx)
  
  def fx(self, x):
    return lognorm.pdf(x, s=self.zetafx, loc=0.00, scale=np.exp(self.lambdafx))
  
  def hx(self, x):
    return lognorm.pdf(x, s=self.zetahx, loc=0.00, scale=np.exp(self.lambdahx))
  
  def zf(self, x):
    return (np.log(x) - self.lambdafx) / self.zetafx