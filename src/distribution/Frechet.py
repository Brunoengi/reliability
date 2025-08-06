from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.base_types.validate_dictionary import ValidateDictionary
from utils.validate.base_types.validate_class import ValidateClass 
import scipy.optimize
import numpy as np
from scipy.special import gamma
import scipy.optimize
import numpy as np
import scipy.linalg
from math import log
from scipy.stats import norm, uniform, lognorm, gumbel_r, invweibull, weibull_min, beta as beta_dist, gamma as gamma_dist, multivariate_normal
from scipy.optimize import fsolve, newton
from scipy.linalg import cholesky
from math import sqrt, pi, log
from scipy.special import gamma

class Frechet(AbstractDistribution):
  def __init__(self, props: dict):
    
    ##self.validate_specific_parameters(props)
    super().__init__(props)
    
    self.deltafx = self.sigmafx / self.mufx
    self.kapa0 = 2.50
    self.gsinal = -1.00
    self.kapaf = scipy.optimize.newton(self.fkapa, self.kapa0, args=(self.deltafx, self.gsinal))
    self.vfn = self.mufx / gamma(1.00 - 1.00 / self.kapaf)
    self.deltahx = self.sigmahx / self.muhx
    self.kapah = scipy.optimize.newton(self.fkapa, self.kapa0, args=(self.deltahx, self.gsinal))
    self.vhn = self.muhx / gamma(1.00 - 1.00 / self.kapah)
    
    ##ValidateClass.has_invalid_key(self, 'varname', 'vardist', 'varmean', 'varcov', 'varhmean','varstd')

  def validate_specific_parameters(self, props):
    ValidateDictionary.is_dictionary(props)
    ValidateDictionary.has_keys(props, 'varmean', 'varcov')
    ValidateDictionary.check_if_exists(props, 'varcov', lambda d, k: ValidateDictionary.is_greater_or_equal_than(d, k, 0))
  
  @staticmethod
  def fkapa(kapa, deltax, gsignal):
        fk = 1.00 + deltax ** 2 - gamma(1.00 + gsignal * 2.00 / kapa) / gamma(1.00 + gsignal * 1.00 / kapa) ** 2
        return fk
  
  def x_uncorrelated(self, ns):
    return invweibull.rvs(c=self.kapah, loc=0.00, scale=self.vhn, size=ns)
  
  def fx_uncorrelated(self, x):
    return invweibull.pdf(x, c=self.kapaf, loc=0.00, scale=self.vfn)
  
  def hx_uncorrelated(self, x):
    return invweibull.pdf(x, c=self.kapah, loc=0.00, scale=self.vhn)
  
  def x_correlated(self, zk_col):
    uk = norm.cdf(zk_col)
    return self.vhn / (np.log(1 / uk)) ** (1 / self.kapah)
  
  def fx_correlated(self, x):
    ynf = x / self.vfn
    return invweibull.pdf(ynf, self.kapaf) / self.vfn
  
  def hx_correlated(self, x):
    ynh = x / self.vhn
    return invweibull.pdf(ynh, self.kapah) / self.vhn
  
  def zf_correlated(self, x):
    ynf = x / self.vfn
    cdfx = invweibull.cdf(ynf, self.kapaf)
    return norm.ppf(cdfx)