from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.base_types.validate_dictionary import ValidateDictionary
from utils.validate.base_types.validate_class import ValidateClass 
import numpy as np
from scipy.optimize import fsolve
from scipy.stats import norm, uniform, lognorm, gumbel_r, invweibull, weibull_min, beta as beta_dist

class Beta(AbstractDistribution):
  def __init__(self, props):

    ##self.validate_specific_parameters(props)

    self.a = props['parameter1']
    self.b = props['parameter2']
    self.q = props['parameter3']
    self.r = props['parameter4']

    self.varmean = float(self.a + self.q / (self.q + self.r) * (self.b - self.a))
    self.varstd = float(np.sqrt((self.q * self.r) / ((self.q + self.r) **2 * (self.q + self.r + 1)) * (self.b - self.a) ** 2))
    props['varmean'] = self.varmean
    props['varstd'] = self.varstd
    
    super().__init__(props)
    
    self.loc = self.a
    self.scale = (self.b - self.a)
    self.ah, self.bh =  fsolve(self.beta_limits, (1, 1), args= ( self.muhx, self.sigmahx, self.q, self.r)) 
    
    self.loch = self.ah
    self.scaleh = (self.bh - self.ah)  
    

    
    ##ValidateClass.has_invalid_key(self, 'varname', 'vardist','varmean','varcov','varstd','varhmean','parameter1','parameter2','parameter3','parameter4')

  def validate_specific_parameters(self, props):
    ValidateDictionary.is_dictionary(props)
    ValidateDictionary.check_possible_arrays_keys(props, ['parameter1', 'parameter2', 'parameter3', 'parameter4'])
    
  def instrumental_properties(self, varhmean):
    self.varhmean = varhmean
    self.muhx = self.varhmean
    self.ah, self.bh = fsolve(self.beta_limits, (1, 1), args=(self.muhx, self.sigmahx, self.q, self.r))
    # Ensure that the instrumental support function covers the original function
    if self.ah > self.loc:
        self.ah = self.loc
    if self.bh < self.loc + self.scale:
        self.bh = self.loc + self.scale
    self.loch = self.ah
    self.scaleh = self.bh - self.ah 
    
  @staticmethod  
  def beta_limits(vars, mux, sigmax, q, r):
    a, b = vars
    eq1 = a + q / (q + r) * (b - a) - mux
    eq2 = ((q * r) / ((q + r) ** 2 * (q + r + 1))) ** (0.50) * (b - a) - sigmax
    return [eq1, eq2]
  
  def x_uncorrelated(self, ns):
    return beta_dist.rvs(self.q, self.r, self.loch, self.scaleh, size=ns)
  
  def fx_uncorrelated(self, x):
    return beta_dist.pdf(x, self.q, self.r, self.loc, self.scale)
  
  def hx_uncorrelated(self, x):
    return beta_dist.pdf(x, self.q, self.r, self.loch, self.scaleh)
  
  def x_correlated(self, zk_col):
    uk = norm.cdf(zk_col)
    return beta_dist.ppf(uk, self.q, self.r, loc=self.loc, scale=self.scale)
  
  def fx_correlated(self, x):
    return beta_dist.pdf(x, self.q, self.r, loc=self.loc, scale=self.scale)
  
  def hx_correlated(self, x):
    return beta_dist.pdf(x, self.q, self.r, loc=self.loch, scale=self.scaleh)
  
  def zf_correlated(self, x):
    cdfx = beta_dist.cdf(x, self.q, self.r, loc=self.loc, scale=self.scale)
    return norm.ppf(cdfx)