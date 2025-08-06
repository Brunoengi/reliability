import scipy.optimize
import numpy as np
import scipy.linalg
from math import log
from scipy.stats import norm, uniform, lognorm, gumbel_r, invweibull, weibull_min, beta as beta_dist, gamma as gamma_dist, multivariate_normal
from scipy.optimize import fsolve, newton
from scipy.linalg import cholesky
from math import sqrt, pi, log
from scipy.special import gamma
from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.base_types.validate_dictionary import ValidateDictionary

class Weibull(AbstractDistribution):
  def __init__(self, props: dict):

    #self.validate_specific_parameters(props)
    super().__init__(props)
    self.varinf = float(props['parameter3'])
    self.epsilon = self.varinf
    self.deltafx = self.sigmafx / (self.mufx - self.epsilon)
    self.kapa0 = 2.50
    self.gsinal = 1.00
    self.kapaf = scipy.optimize.newton(self.fkapa, self.kapa0, args=(self.deltafx, self.gsinal))
    self.w1f = (self.mufx - self.epsilon) / gamma(1.00 + 1.00 / self.kapaf) + self.epsilon
    self.deltahx = self.sigmahx / (self.muhx - self.epsilon)
    self.kapah = scipy.optimize.newton(self.fkapa, self.kapa0, args=(self.deltahx, self.gsinal))
    self.w1h = (self.muhx - self.epsilon) / gamma(1.00 + 1.00 / self.kapah) + self.epsilon
    
    #ValidateClass.has_invalid_key(self, 'varname', 'vardist', 'varmean', 'varstd', 'parameter3', 'varhmean','varinf')

  def validate_specific_parameters(self, props):
    ValidateDictionary.is_dictionary(props)
    ValidateDictionary.has_keys(props, 'varmean', 'varstd', 'parameter3')
    ValidateDictionary.check_if_exists(props, 'varcov', lambda d, k: ValidateDictionary.is_greater_or_equal_than(d, k, 0))
   
  @staticmethod 
  def fkapa(kapa, deltax, gsignal):
        fk = 1.00 + deltax ** 2 - gamma(1.00 + gsignal * 2.00 / kapa) / gamma(1.00 + gsignal * 1.00 / kapa) ** 2
        return fk
      
  def x_uncorrelated(self, ns):
    return weibull_min.rvs(c=self.kapah, loc=self.epsilon, scale=self.w1h-self.epsilon, size=ns)
  
  def fx_uncorrelated(self, x):
    return weibull_min.pdf(x, c=self.kapaf, loc=self.epsilon, scale=self.w1f-self.epsilon)
  
  def hx_uncorrelated(self, x):
    return weibull_min.pdf(x, c=self.kapah, loc=self.epsilon, scale=self.w1h-self.epsilon)
  
  def x_correlated(self, zk_col):
    uk = norm.cdf(zk_col)
    return (self.w1h - self.epsilon) * (np.log(1 / (1 - uk))) ** (1 / self.kapah) + self.epsilon
  
  def fx_correlated(self, x):
    ynf = (x - self.epsilon) / (self.w1f - self.epsilon)
    return weibull_min.pdf(ynf, self.kapaf) / (self.w1f - self.epsilon)
  
  def hx_correlated(self, x):
    ynh = (x - self.epsilon) / (self.w1h - self.epsilon)
    return weibull_min.pdf(ynh, self.kapah) / (self.w1h - self.epsilon)
    
  def zf_correlated(self, x):
    ynf = (x - self.epsilon) / (self.w1f - self.epsilon)
    cdfx = weibull_min.cdf(ynf, self.kapaf)
    return norm.ppf(cdfx)
    