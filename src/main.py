import numpy as np
from scipy.stats import norm
from scipy.stats import lognorm
from scipy.stats import uniform
from scipy.stats import gumbel_r
from scipy.stats import invweibull
from scipy.stats import weibull_min
from scipy.stats import multivariate_normal
from scipy.stats import beta as beta_dist
from scipy.stats import gamma as gamma_dist
import scipy.optimize
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
import scipy.linalg
from scipy.special import gamma
import time
from utils.distribution import createDistribution
from visualize import DataVisualize
from utils.validate.domain_types.validate_dvar import ValidateDvar
from utils.validate.domain_types.validate_gx import ValidateGx
from utils.validate.domain_types.validate_corrmatrix import ValidateCorrelationMatrix
from correlation import Correlation
from methods.transformation import TransformationMethods
from methods.simulation.monte_carlo import MonteCarloMethods


class Reliability():

    def __init__(self, xvar, dvar, gx, x0=None, Rz=None):
        
        #Simple setters
        self.dvar = dvar 
        self.fel = gx
        self.ndvar = len(dvar)
        self.nxvar = len(xvar)

        #Data processing - Complex setters
        self._xvar = self.set_xvar(xvar, x0)
        self._Rz = self.set_Rz(Rz)
        self._d = self.set_d(self.dvar, self.ndvar)
        self._x0 = self.set_x0(x0, self.nxvar, self.xvar)

        self._correlation = Correlation({'Rz': self._Rz, 'xvar': self._xvar, 'nxvar': self.nxvar})
        self._transformation =  TransformationMethods(self)
        self._simulation = MonteCarloMethods(self)
           
    @property
    def xvar(self):
      return self._xvar

    @property
    def dvar(self):
      return self._dvar

    @property
    def gx(self):
      return self._gx

    @property
    def Rz(self):
      return self._Rz

    @property
    def x0(self):
      return self._x0

    @property
    def nxvar(self):
      return len(self.xvar)
    
    @property
    def ndvar(self):
      return len(self.dvar)
    
    @property
    def d(self):
        return self._d
    
    @property
    def correlation(self):
      return self._correlation

    @xvar.setter
    def xvar(self, xvar):
        self._xvar = xvar

    @dvar.setter
    def dvar(self, dvar): 
      ValidateDvar(dvar)
      self._dvar = dvar

    @gx.setter
    def gx(self, gx): 
      ValidateGx(gx)
      self._gx = gx

    @Rz.setter
    def Rz(self, Rz):
        self._Rz = Rz

    @x0.setter
    def x0(self, x0):
        self._x0 = x0 

    def set_xvar(self, xvar, x0): 
      nxvar = len(xvar)
      for var in xvar:              
        #
        # Setting standard variable distribution names, create distribution and update the `var` dictionary with all distribution attributes
        #
        var.update(vars(createDistribution(var)))
        #
        # Initial values of the aleatory variables
        #
        if x0 is None:
            #
            # Original mean of the variables x
            #
            i = -1
            x0 = np.zeros(nxvar)
            for var in xvar:
                i += 1
                # Mean value of the random variables x
                x0[i] = float(var['varmean'])
                var['varhmean'] = float(var['varmean'])
        else:
            i = -1
            for var in xvar:
                i += 1
                # Mean value of the random variables x
                var['varhmean'] = x0[i]

      return xvar

    def set_Rz(self, Rz):
      if Rz is None:
          return np.eye(self.nxvar)
      else:
        ValidateCorrelationMatrix(Rz)
        return np.array(Rz)
        
    @nxvar.setter
    def nxvar(self, nxvar): 
      self._nxvar = nxvar
  
    @ndvar.setter
    def ndvar(self, ndvar): 
      self._ndvar = ndvar

    @correlation.setter
    def correlation(self, correlation): 
      self._correlation = correlation

    def set_d(self, dvar, ndvar):
      #
      # Initial values of the design variables
      #
      i = -1
      d = np.zeros(ndvar)
      for var in dvar:
          i += 1
          # Mean value of the random variables x
          d[i] = float(var['varvalue'])
      
      #
      # Setting variables initial values
      #
      return d

    def set_x0(self, x0, nxvar, xvar):
      if x0 is None:
          #
          # Original mean of the variables x
          #
          i = -1
          x0 = np.zeros(nxvar)
          for var in xvar:
              i += 1
              # Mean value of the random variables x
              x0[i] = float(var['varmean'])
              var['varhmean'] = float(var['varmean'])
      else:
          i = -1
          for var in xvar:
              i += 1
              # Mean value of the random variables x
              var['varhmean'] = x0[i]

      return x0

    def __getattr__(self, name):
      for module in [self._transformation, self._simulation]:
        if hasattr(module, name):
          return getattr(module, name)
      raise AttributeError(f"'Reliability' object has no attribute '{name}'")

    def var_gen(self, ns, indexes_correlated_xvar, nsigma=1.00, iprint=False):
        """

           Random variables generator for the Monte Carlo Simulation methods

        """

        #Get only correlated variables
        xvar_correlated = [self.xvar[i] for i in indexes_correlated_xvar]
        nxvar_correlated = len(xvar_correlated)


        # Generation of uniform random numbers
        uk_cycle = np.random.rand(ns, nxvar_correlated)
        #


        def fkapa(kapa, deltax, gsignal):
            fk = 1.00 + deltax ** 2 - gamma(1.00 + gsignal * 2.00 / kapa) / gamma(1.00 + gsignal * 1.00 / kapa) ** 2
            return fk
        
        def beta_limits(vars, mux, sigmax, q, r):
            a, b = vars
            eq1 = a + q / (q + r) * (b - a) - mux
            eq2 = ((q * r) / ((q + r) ** 2 * (q + r + 1))) ** (0.50) * (b - a) - sigmax
            return [eq1, eq2]
        
        def uniform_limits(vars, mux, sigmax):
            a, b = vars
            eq1 = (a + b) / 2 - mux
            eq2 = (b - a) / np.sqrt(12.) - sigmax
            return [eq1, eq2]


        x = np.zeros((ns, nxvar_correlated))
        weight = np.ones(ns)
        fx = np.zeros(ns)
        hx = np.zeros(ns)
        fxixj = np.ones(ns)
        yk = np.zeros((ns, nxvar_correlated))



        # Get a sub-matrix only correlated variables
        matrix = self.correlation.Rz_rectify[np.ix_(indexes_correlated_xvar, indexes_correlated_xvar)]


        # Step 1 - Determination of equivalent correlation coefficients and
        #          Jacobian matrix Jzy
        #
        #
        # Cholesky decomposition of the correlation matrix
        #
        L = scipy.linalg.cholesky(matrix, lower=True)
        Jzy = np.copy(L)

        #
        # Generation of Gaussian correlated random numbers
        #

        yk = norm.ppf(uk_cycle)
        zf = np.zeros((ns, nxvar_correlated))
        zk = np.dot(Jzy, yk.T).T
        # zk = multivariate_normal.rvs(mean = np.zeros(self.nxvar), cov = self.Rz, size=ns)


        i = -1
        for var in xvar_correlated:
            i += 1
            
            if var['varstd'] == 0.00:
                var['varstd'] = float(var['varcov']) * float(var['varmean'])
            #
            #
            # Normal distribution
            #
            namedist = var['vardist']
            if namedist.lower() == 'gauss':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                x[:, i] = muhx + sigmahx * zk[:, i]
                fx = norm.pdf(x[:, i], mufx, sigmafx)
                hx = norm.pdf(x[:, i], muhx, sigmahx)
                zf[:, i] = (x[:, i]-mufx)/sigmafx
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)
            #
            # Uniform or constant distribution
            #
            
            elif namedist.lower() == 'uniform':
                a = float(var['parameter1'])
                b = float(var['parameter2'])
                
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                ah, bh =  fsolve(uniform_limits, (1, 1), args= (muhx, sigmahx))  
                                
                uk = norm.cdf(zk[:, i])
                x[:, i] = ah + (bh - ah) * uk
                zf[:, i] = norm.ppf(uk)
                fx = uniform.pdf(x[:, i], a, b-a)
                hx = uniform.pdf(x[:, i], ah, bh-ah)
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)
            #
            # Lognormal distribution
            #
            elif namedist.lower() == 'lognorm':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                zetafx = np.sqrt(np.log(1.00 + (sigmafx / mufx) ** 2))
                lambdafx = np.log(mufx) - 0.5 * zetafx ** 2
                zetahx = np.sqrt(np.log(1.00 + (sigmahx / muhx) ** 2))
                lambdahx = np.log(muhx) - 0.5 * zetahx ** 2
                x[:, i] = np.exp(lambdahx + zk[:, i] * zetahx)
                zf[:, i] = (np.log(x[:, i])-lambdafx) / zetafx
                fx = norm.pdf(np.log(x[:, i]), lambdafx, zetafx)
                hx = norm.pdf(np.log(x[:, i]), lambdahx, zetahx)
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)

            #
            # Gumbel distribution
            #
            elif namedist.lower() == 'gumbel':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                alphafn = np.pi / np.sqrt(6.00) / sigmafx
                ufn = mufx - np.euler_gamma / alphafn
                betafn = 1.00 / alphafn
                alphahn = np.pi / np.sqrt(6.00) / sigmahx
                uhn = muhx - np.euler_gamma / alphahn
                betahn = 1.00 / alphahn
                uk = norm.cdf(zk[:, i])
                x[:, i] = uhn - betahn * np.log(np.log(1. / uk))
                cdfx = gumbel_r.cdf(x[:, i], ufn, betafn)
                zf[:, i] = norm.ppf(cdfx, 0, 1)
                fx = gumbel_r.pdf(x[:, i], ufn, betafn)
                hx = gumbel_r.pdf(x[:, i], uhn, betahn)
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)

            #
            # Frechet distribution
            #
            elif namedist.lower() == 'frechet':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                deltafx = sigmafx / mufx
                kapa0 = 2.50
                gsinal = -1.00
                kapaf = scipy.optimize.newton(fkapa, kapa0, args=(deltafx, gsinal))
                vfn = mufx / gamma(1.00 - 1.00 / kapaf)
                deltahx = sigmahx / muhx
                kapa0 = 2.50
                gsinal = -1.00
                kapah = scipy.optimize.newton(fkapa, kapa0, args=(deltahx, gsinal))
                vhn = muhx / gamma(1.00 - 1.00 / kapah)
                uk = norm.cdf(zk[:, i])
                x[:, i] = vhn / (np.log(1. / uk)) ** (1. / kapah)
                ynf = x[:, i] / vfn
                ynh = x[:, i] / vhn
                cdfx = invweibull.cdf(ynf, kapaf)
                zf[:, i] = norm.ppf(cdfx, 0, 1)
                fx = invweibull.pdf(ynf, kapaf) / vfn
                hx = invweibull.pdf(ynh, kapah) / vhn
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)

            #
            #
            # Weibull distribution - minimum
            #
            elif namedist.lower() == 'weibull':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                epsilon = float(var['varinf'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                deltafx = sigmafx / (mufx - epsilon)
                kapa0 = 2.50
                gsinal = 1.00
                kapaf = scipy.optimize.newton(fkapa, kapa0, args=(deltafx, gsinal))
                w1f = (mufx - epsilon) / gamma(1.00 + 1.00 / kapaf) + epsilon
                deltahx = sigmahx / (muhx - epsilon)
                kapa0 = 2.50
                gsinal = 1.00
                kapah = scipy.optimize.newton(fkapa, kapa0, args=(deltahx, gsinal))
                w1h = (muhx - epsilon) / gamma(1.00 + 1.00 / kapah) + epsilon
                uk = norm.cdf(zk[:, i])
                x[:, i] = (w1h - epsilon) * (np.log(1./(1. - uk))) ** (1. / kapah) + epsilon
                ynf = (x[:, i] - epsilon) / (w1f - epsilon)
                ynh = (x[:, i] - epsilon) / (w1h - epsilon)
                cdfx = weibull_min.cdf(ynf, kapaf)
                zf[:, i] = norm.ppf(cdfx, 0, 1)
                fx = weibull_min.pdf(ynf, kapaf) / (w1f - epsilon)
                hx = weibull_min.pdf(ynh, kapah) / (w1h - epsilon)
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)

            #
            #
            # Beta distribution
            #
            elif namedist.lower() == 'beta':
                a = float(var['parameter1'])
                b = float(var['parameter2'])
                q = float(var['parameter3'])
                r = float(var['parameter4'])
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                loc = a
                scale = (b - a)
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                ah, bh =  fsolve(beta_limits, (1, 1), args= ( muhx, sigmahx, q, r))  
                loch = ah
                scaleh = (bh - ah)        
                uk = norm.cdf(zk[:, i])
                x[:, i] = beta_dist.ppf(uk, q, r, loc, scale)
                fx = beta_dist.pdf(x[:, i], q, r, loc, scale)
                hx = beta_dist.pdf(x[:, i], q, r, loch, scaleh)
                cdfx = beta_dist.cdf(x[:, i], q, r, loc, scale)
                zf[:, i] = norm.ppf(cdfx, 0, 1)
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)

            #
            #
            # Gamma distribution
            #
            elif namedist.lower() == 'gamma':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                deltafx = sigmafx / mufx
                k = 1. / deltafx ** 2
                v = k / mufx
                a = k
                loc = 0.00
                scale = 1. / v
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                deltahx = sigmahx / muhx
                kh = 1. / deltahx ** 2
                vh = kh / muhx
                ah = kh
                loch = 0.00
                scaleh = 1. / vh
                uk = norm.cdf(zk[:, i])
                x[:, i] = gamma_dist.ppf(uk, ah, loch, scaleh)
                fx = gamma_dist.pdf(x[:, i], a, loc, scale)
                hx = gamma_dist.pdf(x[:, i], ah, loch, scaleh)
                cdfx = gamma_dist.cdf(x[:, i], a, loc, scale)
                zf[:, i] = norm.ppf(cdfx, 0, 1)
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)
                
        
        norm_multivarf = multivariate_normal(mean=None, cov=matrix)
        phif = list(map(norm_multivarf.pdf, zf))
        phif = np.array(phif)
        norm_multivarh = multivariate_normal(mean=None, cov=matrix)
        phih = list(map(norm_multivarh.pdf, zk))
        phih = np.array(phih)
        weight = weight * phif / phih
        fxixj = fxixj * phif

        return x, weight, fxixj
    
    def var_rvs(self, ns, indexes_uncorrelated_xvar, nsigma=1.00, iprint=False):
        """

           Random variables generator for the Monte Carlo Simulation methods

        """

        #Get only correlated variables
        xvar_uncorrelated = [self.xvar[i] for i in indexes_uncorrelated_xvar]
        nxvar_uncorrelated = len(xvar_uncorrelated)

        def fkapa(kapa, deltax, gsignal):
            fk = 1.00 + deltax ** 2 - gamma(1.00 + gsignal * 2.00 / kapa) / gamma(1.00 + gsignal * 1.00 / kapa) ** 2
            return fk
        
        def beta_limits(vars, mux, sigmax, q, r):
            a, b = vars
            eq1 = a + q / (q + r) * (b - a) - mux
            eq2 = ((q * r) / ((q + r) ** 2 * (q + r + 1))) ** (0.50) * (b - a) - sigmax
            return [eq1, eq2]
        
        def uniform_limits(vars, mux, sigmax):
            a, b = vars
            eq1 = (a + b) / 2 - mux
            eq2 = (b - a) / np.sqrt(12.) - sigmax
            return [eq1, eq2]

        

        x = np.zeros((ns, nxvar_uncorrelated))
        weight = np.ones(ns)
        fx = np.zeros(ns)
        hx = np.zeros(ns)
        fxixj = np.ones(ns)
        
        #
        # Step 1 - Determination of equivalent correlation coefficients and
        #          Jacobian matrix Jzy
        #
        #
        # Cholesky decomposition of the correlation matrix
        #

        L = scipy.linalg.cholesky(self.correlation.Rz_rectify, lower=True)
        Jzy = np.copy(L)

        #
        # Generation of Gaussian random numbers
        #

        
        zf = np.zeros((ns, nxvar_uncorrelated))
        
        #
        i = -1
        for var in xvar_uncorrelated:
            i += 1
            if var['varstd'] == 0.00:
                var['varstd'] = float(var['varcov']) * float(var['varmean'])
            if iprint:
                print(self.xvar[i])
            #
            #
            # Normal distribution
            #
            namedist = var['vardist']
            if namedist.lower() == 'gauss':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                x[:, i] = norm.rvs(loc=muhx, scale=sigmahx, size=ns)
                fx = norm.pdf(x[:, i], mufx, sigmafx)
                hx = norm.pdf(x[:, i], muhx, sigmahx)
                weight = weight * (fx / hx)
                fxixj = fxixj * fx 
            #
            # Uniform or constant distribution
            #
            
            elif namedist.lower() == 'uniform':
                a = float(var['parameter1'])
                b = float(var['parameter2'])
                
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                ah, bh =  fsolve(uniform_limits, (1, 1), args= (muhx, sigmahx))  
                                
                
                x[:, i] = uniform.rvs(loc=ah, scale= (bh-ah), size = ns)
                fx = uniform.pdf(x[:, i], a, b-a)
                hx = uniform.pdf(x[:, i], ah, bh-ah)
                weight = weight * (fx / hx)
                fxixj = fxixj * fx 
            #
            # Lognormal distribution
            #
            elif namedist.lower() == 'lognorm':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                zetafx = np.sqrt(np.log(1.00 + (sigmafx / mufx) ** 2))
                lambdafx = np.log(mufx) - 0.5 * zetafx ** 2
                zetahx = np.sqrt(np.log(1.00 + (sigmahx / muhx) ** 2))
                lambdahx = np.log(muhx) - 0.5 * zetahx ** 2
                x[:, i] = lognorm.rvs(s=zetahx, loc=0.00, scale=np.exp(lambdahx), size=ns)
                fx = lognorm.pdf(x[:, i], s=zetafx, loc=0.00, scale=np.exp(lambdafx))
                hx = lognorm.pdf(x[:, i], s=zetahx, loc=0.00, scale=np.exp(lambdahx))
                weight = weight * (fx / hx)
                fxixj = fxixj * fx 

            #
            # Gumbel distribution
            #
            elif namedist.lower() == 'gumbel':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                alphafn = np.pi / np.sqrt(6.00) / sigmafx
                ufn = mufx - np.euler_gamma / alphafn
                betafn = 1.00 / alphafn
                alphahn = np.pi / np.sqrt(6.00) / sigmahx
                uhn = muhx - np.euler_gamma / alphahn
                betahn = 1.00 / alphahn
                x[:, i] = gumbel_r.rvs( loc=uhn, scale=betahn, size=ns)
                fx = gumbel_r.pdf(x[:, i], ufn, betafn)
                hx = gumbel_r.pdf(x[:, i], uhn, betahn)
                weight = weight * (fx / hx)
                fxixj = fxixj * fx 

            #
            # Frechet distribution
            #
            elif namedist.lower() == 'frechet':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                deltafx = sigmafx / mufx
                kapa0 = 2.50
                gsinal = -1.00
                kapaf = scipy.optimize.newton(fkapa, kapa0, args=(deltafx, gsinal))
                vfn = mufx / gamma(1.00 - 1.00 / kapaf)
                deltahx = sigmahx / muhx
                kapa0 = 2.50
                gsinal = -1.00
                kapah = scipy.optimize.newton(fkapa, kapa0, args=(deltahx, gsinal))
                vhn = muhx / gamma(1.00 - 1.00 / kapah)
                x[:, i] = invweibull.rvs(c=kapah, loc=0.00, scale=vhn, size=ns)
                fx = invweibull.pdf(x[:, i], c=kapaf, loc=0.00, scale=vfn)
                hx = invweibull.pdf(x[:, i], c=kapah, loc=0.00, scale=vhn)
                weight = weight * (fx / hx)
                fxixj = fxixj * fx 

            #
            #
            # Weibull distribution - minimum
            #
            elif namedist.lower() == 'weibull':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                epsilon = float(var['varinf'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                deltafx = sigmafx / (mufx - epsilon)
                kapa0 = 2.50
                gsinal = 1.00
                kapaf = scipy.optimize.newton(fkapa, kapa0, args=(deltafx, gsinal))
                w1f = (mufx - epsilon) / gamma(1.00 + 1.00 / kapaf) + epsilon
                deltahx = sigmahx / (muhx - epsilon)
                kapa0 = 2.50
                gsinal = 1.00
                kapah = scipy.optimize.newton(fkapa, kapa0, args=(deltahx, gsinal))
                w1h = (muhx - epsilon) / gamma(1.00 + 1.00 / kapah) + epsilon
                x[:, i] = weibull_min.rvs(c=kapah, loc=epsilon, scale=w1h-epsilon, size=ns)
                fx = weibull_min.pdf(x[:, i], c=kapaf, loc=epsilon, scale=w1f-epsilon)
                hx = weibull_min.pdf(x[:, i], c=kapah, loc=epsilon, scale=w1h-epsilon)
                weight = weight * (fx / hx)
                fxixj = fxixj * fx 

            #
            #
            # Beta distribution
            #
            elif namedist.lower() == 'beta':
                a = float(var['parameter1'])
                b = float(var['parameter2'])
                q = float(var['parameter3'])
                r = float(var['parameter4'])
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                loc = a
                scale = (b - a)
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                ah, bh =  fsolve(beta_limits, (1, 1), args= ( muhx, sigmahx, q, r))  
                loch = ah
                scaleh = (bh - ah)        
                x[:, i] = beta_dist.rvs(q, r, loch, scaleh, size=ns)
                fx = beta_dist.pdf(x[:, i], q, r, loc, scale)
                hx = beta_dist.pdf(x[:, i], q, r, loch, scaleh)
                weight = weight * (fx / hx)
                fxixj = fxixj * fx 

            #
            #
            # Gamma distribution
            #
            elif namedist.lower() == 'gamma':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                deltafx = sigmafx / mufx
                k = 1. / deltafx ** 2
                v = k / mufx
                a = k
                loc = 0.00
                scale = 1. / v
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                deltahx = sigmahx / muhx
                kh = 1. / deltahx ** 2
                vh = kh / muhx
                ah = kh
                loch = 0.00
                scaleh = 1. / vh
                x[:, i] = gamma_dist.rvs(ah, loch, scaleh, size=ns)
                fx = gamma_dist.pdf(x[:, i], a, loc, scale)
                hx = gamma_dist.pdf(x[:, i], ah, loch, scaleh)
                weight = weight * (fx / hx)
                fxixj = fxixj * fx 
                

        return x, weight, fxixj

    def multig(self, xvar, dvar, glist, iprint=True):
        """
        Solution of the problem of the reliability of serial system with multiple limit state functions
        According to:
        BECK, A.T.
        Confiabilidade e Segurança das Estruturas
        Rio de Janeiro, Elsevier, 2019.
        """
        ng = int(len(glist))
        nvar = int(len(xvar))
        beta = np.zeros(ng)
        alpha = np.zeros((ng, nvar))
        pf = np.zeros(ng)
        pa = np.zeros((ng, ng))
        pb = np.zeros((ng, ng))
        pfij_inf = np.zeros((ng, ng))
        pfij_sup = np.zeros((ng, ng))

        i = -1
        for gfunction in glist:
            i += 1
            #
            # FORM method for  the multiple g(x) functions
            #
            test = Reliability(xvar, dvar, gfunction, None, None)
            beta[i], x0, alpha[i, :], normgradyk, sigmaxneqk = test.form(iHLRF=True, toler=1.e-3)
            pf[i] = norm.cdf(-beta[i])
        #

        pfinf = pf.max()
        pfsup = pf.sum()
        alpha_sign = np.sign(alpha)
        alpha2 = alpha_sign * alpha ** 2

        #
        # Sort arrays pf, beta, alpha, alpha2, glist in decrescent order of probability of failure
        #
        glist = np.array(glist)
        ig = (-pf).argsort()
        pf = pf[ig]
        beta = beta[ig]
        alpha = alpha[ig, :]
        alpha2 = alpha2[ig, :]
        glist = glist[ig]


        #
        # Calculation of the correlation coefficients between the limit state functions
        #

        ro = np.dot(alpha, alpha.T)

        #
        # Calculation of the p(Aij) and p(Bij) matrices
        #

        for i in range(ng):
            for j in range(ng):
                if i != j:
                    pa[i, j] = norm.cdf(-beta[i]) * norm.cdf(
                        -((beta[j] - ro[i, j] * beta[i]) / np.sqrt(1. - ro[i, j] ** 2)))
                    pb[i, j] = norm.cdf(-beta[j]) * norm.cdf(
                        -((beta[i] - ro[i, j] * beta[j]) / np.sqrt(1. - ro[i, j] ** 2)))


        #
        # Calculation of the terms P(Fi.Fj)sup and P(Fi.Fj)inf
        #

        for i in range(ng):
            for j in range(ng):
                if i != j:
                    if ro[i, j] >= 0.00:
                        pfij_inf[i, j] = pa[i, j] + pb[i, j]
                        pfij_sup[i, j] = np.max([pa[i, j], pb[i, j]])
                    else:
                        pfij_inf[i, j] = np.min([pa[i, j], pb[i, j]])
                        pfij_sup[i, j] = 0.00

        #
        # Calculation of inferior and superior limits for the probability of failure of the system
        #
        #
        # Inferior limit: pf_inf
        #
        pf_inf = pf[0]
        for i in range(1, ng, 1):
            pfi_inf = pf[i]
            for j in range(0, i, 1):
                pfi_inf -= pfij_inf[i, j]
            pf_inf += np.max([0, pfi_inf])

        #
        # Superior limit: pf_sup
        #
        pf_sup = sum(pf)
        for i in range(1, ng, 1):
            pf_sup -= np.max(pfij_sup[i, 0:i], axis=0)

        beta_sup = -norm.ppf(pf_inf)
        beta_inf = -norm.ppf(pf_sup)
        glist = list(glist)

        if iprint:
            #
            # Print the initial results
            #
            print('Initial results:')
            print('pf =', pf)
            print('beta =', beta)
            print('pfinf =', pfinf)
            print('pfsup =', pfsup)
            print('alpha =', alpha)
            print(('alpha2 =', alpha2))
            print('ro =', ro)
            print('pa =', pa)
            print('pb =', pb)
            #
            # Print final results
            #
            print('Final results:')
            print('g list =', glist)
            print('pf =', pf)
            print('beta =', beta)
            print('alpha =', alpha)
            print(('alpha2 =', alpha2))
            print('pf_inf =', pf_inf)
            print('pf_sup =', pf_sup)
            print('beta_inf =', beta_inf)
            print('beta_sup =', beta_sup)

    def generator(self, ns, nsigma=1.00, iprint=False):
            """
            Method to generate random variables

            """
            # Number of variables of the problem
            ns = int(ns)
            
            # Get index of correlated and uncorrelated variables
            index_correlated, index_uncorrelated = self.correlation.correlation_summary()
            total_index = len(index_correlated) + len(index_uncorrelated) 

            print('correlated', index_correlated)
            print('uncorrelated', index_uncorrelated)

            # Empty matrix 
            x = np.empty((ns, total_index))

            #
            # Standard deviation multiplier for MC-IS
            
            #
            # Step 1 - Generation of the random numbers according to their appropriate distribution
            #

            wpc = np.ones(ns)
            fxc = np.ones(ns)
            wpu = np.ones(ns)
            fxu = np.ones(ns)
            
            if index_correlated:
              xpc, wpc, fxc = self.var_gen(ns, index_correlated, nsigma)
            
            if index_uncorrelated:
              xpu, wpu, fxu = self.var_rvs(ns, index_uncorrelated, nsigma)
            


            for i, idx in enumerate(index_correlated):
              x[:, idx] = xpc[:, i]

            # Preencher as colunas dos não correlacionados
            if index_correlated:
              for i, idx in enumerate(index_correlated):
                x[:, idx] = xpc[:, i]

            if index_uncorrelated:
              for i, idx in enumerate(index_uncorrelated):
                x[:, idx] = xpu[:, i]

            wp = wpc * wpu
            print('wp',wp)
          
            fx = fxc * fxu
            print('fx',fx)
            return x, wp, fx
    
             