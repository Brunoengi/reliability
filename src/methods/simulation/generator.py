import scipy.optimize
import numpy as np
import scipy.linalg
from math import log
from scipy.stats import norm, uniform, lognorm, gumbel_r, invweibull, weibull_min, beta as beta_dist, gamma as gamma_dist, multivariate_normal
from scipy.optimize import fsolve, newton
from scipy.linalg import cholesky
from math import sqrt, pi, log
from scipy.special import gamma


class RandomVariablesGenerator:
  def __init__(self, parent):
      ## Get all properties about Monte Carlo Methods
      self.parent = parent

      ## Get all properties about Reliability
      self.reliability = parent.reliability

  def main(self, ns, nsigma=1.00, iprint=False):
        """
        Method to generate random variables, check correlation matrix to use var_gen (correlated variables) or var_rvs (uncorrelated variables)

        """
        # Number of variables of the problem
        ns = int(ns)
        
        # Get index of correlated and uncorrelated variables
        index_correlated, index_uncorrelated = self.reliability.correlation.correlation_summary()
        total_index = len(index_correlated) + len(index_uncorrelated) 

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
        fx = fxc * fxu

        return x, wp, fx

  def var_gen(self, ns, indexes_correlated_xvar, nsigma=1.00, iprint=False):
      """
      Random variables generator for Monte Carlo Simulation methods, only to correlated variables
      """

      xvar_correlated = [self.reliability.xvar[i] for i in indexes_correlated_xvar]
      nxvar_correlated = len(xvar_correlated)

      # Generate initial matrices and variables
      uk_cycle = np.random.rand(ns, nxvar_correlated)
      x = np.zeros((ns, nxvar_correlated))
      weight = np.ones(ns)
      fxixj = np.ones(ns)
      zf = np.zeros((ns, nxvar_correlated))

      # Correlation submatrix for correlated variables
      matrix = self.reliability.correlation.Rz_rectify[np.ix_(indexes_correlated_xvar, indexes_correlated_xvar)]

      # Cholesky to generate correlated Gaussian samples
      L = scipy.linalg.cholesky(matrix, lower=True)
      yk = norm.ppf(uk_cycle)       # Independent standard normals
      zk = yk @ L.T                 # Correlated normals

      def update_weights(fx, hx, zf_col, zk_col):
          # Update weights and fxixj product for each variable
          standard_norm_pdf_zf = norm.pdf(zf_col)
          standard_norm_pdf_zk = norm.pdf(zk_col)
          w = (fx / standard_norm_pdf_zf) / (hx / standard_norm_pdf_zk)
          
          print('standard_norm_pdf_zf',standard_norm_pdf_zf)
          print('standard_norm_pdf_zk', standard_norm_pdf_zk)
          print('fx', fx)
          print('hx', hx)
          print('w', w)
          
          return w, fx / standard_norm_pdf_zf

      for i, var in enumerate(xvar_correlated):
          # Adjust std if needed
          if var['varstd'] == 0.0:
              var['varstd'] = float(var['varcov']) * float(var['varmean'])

          namedist = var['vardist'].lower()
          mufx = float(var['varmean'])
          sigmafx = float(var['varstd'])
          muhx = float(var['varhmean'])
          sigmahx = nsigma * sigmafx
          
          ##Tem que arrumar no método de Bucher e Adaptive
          print('mufx', float(var['varmean']), self.reliability.xvarClassCorrelated[i].mufx)
          print('sigmafx', float(var['varstd']), self.reliability.xvarClassCorrelated[i].sigmafx)
          print('muhx', float(var['varhmean']), self.reliability.xvarClassCorrelated[i].muhx)
          print('sigmahx', nsigma * sigmafx, self.reliability.xvarClassCorrelated[i].sigmahx)

          zk_col = zk[:, i]

          if namedist == 'gauss':
              x[:, i] = self.reliability.xvarClassCorrelated[i].x_correlated(zk_col)
              fx = self.reliability.xvarClassCorrelated[i].fx(x[:, i])
              hx = self.reliability.xvarClassCorrelated[i].hx(x[:, i]) 
              zf[:, i] = self.reliability.xvarClassCorrelated[i].zf(x[:, i])

          elif namedist == 'uniform':
              x[:, i] = self.reliability.xvarClassCorrelated[i].x_correlated(zk_col)
              zf[:, i] = self.reliability.xvarClassCorrelated[i].zf()
              fx = self.reliability.xvarClassCorrelated[i].fx(x[:, i])
              hx = self.reliability.xvarClassCorrelated[i].hx(x[:, i]) 

          elif namedist == 'lognorm':
            #distribution = self.reliability.xvarClass[i]
            
            
            # print('mufx     ', mufx,      distribution.mufx)
            # print('sigmafx  ', sigmafx,   distribution.sigmafx)
            # print('muhx     ', muhx,      distribution.muhx)
            # print('sigmahx  ', sigmahx,   distribution.sigmahx)
            # print('zetafx   ', zetafx,    distribution.zetafx)
            # print('lambdafx ', lambdafx,  distribution.lambdafx)
            # print('zetahx   ', zetahx,    distribution.zetahx)
            # print('lambdahx ', lambdahx,  distribution.lambdahx)
                    
            
            x[:, i] = self.reliability.xvarClassCorrelated[i].x_correlated(zk_col)
            fx = self.reliability.xvarClassCorrelated[i].fx_correlated(x[:, i])
            hx = self.reliability.xvarClassCorrelated[i].hx_correlated(x[:, i])
            zf[:, i] = self.reliability.xvarClassCorrelated[i].zf_correlated(x[:, i])
            

          elif namedist == 'gumbel':
              
            #   alphafn = pi / sqrt(6) / sigmafx
            #   euler_gamma = 0.5772156649015329
            #   ufn = mufx - euler_gamma / alphafn
            #   betafn = 1 / alphafn
            #   alphahn = pi / sqrt(6) / sigmahx
            #   uhn = muhx - euler_gamma / alphahn
            #   betahn = 1 / alphahn
            #   print('alphafn', alphafn, self.reliability.xvarClass[i].alphafn)
            #   print('ufn', ufn, self.reliability.xvarClass[i].ufn)
            #   print('betafn', betafn, self.reliability.xvarClass[i].betafn)
            #   print('alphahn', alphahn, self.reliability.xvarClass[i].alphahn)
            #   print('uhn', uhn, self.reliability.xvarClass[i].uhn)
            #   print('betahn', betahn, self.reliability.xvarClass[i].betahn)
              
              
              uk = norm.cdf(zk_col)
              x[:, i] = self.reliability.xvarClassCorrelated[i].x_correlated(zk_col)
              zf[:, i] = self.reliability.xvarClassCorrelated[i].zf(zk_col)
              fx = self.reliability.xvarClassCorrelated[i].fx(x[:, i])
              hx = self.reliability.xvarClassCorrelated[i].hx(x[:, i])

          elif namedist == 'frechet':
            #   deltafx = sigmafx / mufx
              def fkapa(kapa, delta, gsignal):
                  return 1 + delta**2 - gamma(1 + 2 * gsignal / kapa) / (gamma(1 + gsignal / kapa) ** 2)

            #   kapaf = newton(fkapa, 2.5, args=(deltafx, -1))
            #   vfn = mufx / gamma(1 - 1 / kapaf)
            #   deltahx = sigmahx / muhx
            #   kapah = newton(fkapa, 2.5, args=(deltahx, -1))
            #   vhn = muhx / gamma(1 - 1 / kapah)
              
              
              
              x[:, i] = self.reliability.xvarClassCorrelated[i].x_correlated(zk_col)
              zf[:, i] = self.reliability.xvarClassCorrelated[i].zf_correlated(x[:, i])
              fx = self.reliability.xvarClassCorrelated[i].fx_correlated(x[:, i])
              hx = self.reliability.xvarClassCorrelated[i].hx_correlated(x[:, i])

          elif namedist == 'weibull':
            #   def fkapa(kapa, delta, gsignal):
            #       return 1 + delta**2 - gamma(1 + 2 * gsignal / kapa) / (gamma(1 + gsignal / kapa) ** 2)
            #   epsilon = float(var['varinf'])
            #   deltafx = sigmafx / (mufx - epsilon)
            #   kapaf = newton(fkapa, 2.5, args=(deltafx, 1))
            #   w1f = (mufx - epsilon) / gamma(1 + 1 / kapaf) + epsilon
            #   deltahx = sigmahx / (muhx - epsilon)
            #   kapah = newton(fkapa, 2.5, args=(deltahx, 1))
            #   w1h = (muhx - epsilon) / gamma(1 + 1 / kapah) + epsilon
              
            #   uk = norm.cdf(zk_col)
              
              x[:, i] = self.reliability.xvarClassCorrelated[i].x_correlated(zk_col)
            #   ynf = (x[:, i] - epsilon) / (w1f - epsilon)
            #   ynh = (x[:, i] - epsilon) / (w1h - epsilon)
            #   cdfx = weibull_min.cdf(ynf, kapaf)
              zf[:, i] = self.reliability.xvarClassCorrelated[i].zf_correlated(x[:, i])
              fx = self.reliability.xvarClassCorrelated[i].fx_correlated(x[:, i])
              hx = self.reliability.xvarClassCorrelated[i].hx_correlated(x[:, i])

          elif namedist == 'beta':
            #   a = float(var['parameter1'])
            #   b = float(var['parameter2'])
            #   q = float(var['parameter3'])
            #   r = float(var['parameter4'])

            #   def beta_limits(vars, mux, sigmax, q, r):
            #       a_, b_ = vars
            #       eq1 = a_ + q / (q + r) * (b_ - a_) - mux
            #       eq2 = sqrt((q * r) / (((q + r) ** 2) * (q + r + 1))) * (b_ - a_) - sigmax
            #       return [eq1, eq2]

            #   ah, bh = fsolve(beta_limits, (a, b), args=(muhx, sigmahx, q, r))
              
            #   loc, scale = a, b - a
            #   loch, scaleh = ah, bh - ah
              
            #   uk = norm.cdf(zk_col)
              x[:, i] = self.reliability.xvarClassCorrelated[i].x_correlated(zk_col) 
              fx = self.reliability.xvarClassCorrelated[i].fx_correlated(x[:, i])
              hx = self.reliability.xvarClassCorrelated[i].hx_correlated(x[:, i])
              zf[:, i] = self.reliability.xvarClassCorrelated[i].zf_correlated(x[:, i])

          elif namedist == 'gamma':
            #   deltafx = sigmafx / mufx
              
              
            #   k = 1 / deltafx ** 2
            #   v = k / mufx
            #   a, loc, scale = k, 0, 1 / v
            #   deltahx = sigmahx / muhx
            #   kh = 1 / deltahx ** 2
            #   vh = kh / muhx
            #   ah, loch, scaleh = kh, 0, 1 / vh

            #   uk = norm.cdf(zk_col)
              x[:, i] = self.reliability.xvarClassCorrelated[i].x_correlated(zk_col) 
              fx = self.reliability.xvarClassCorrelated[i].fx_correlated(x[:, i])
              hx = self.reliability.xvarClassCorrelated[i].hx_correlated(x[:, i])
              
            #   cdfx = gamma_dist.cdf(x[:, i], a, loc=loc, scale=scale)
              zf[:, i] = self.reliability.xvarClassCorrelated[i].zf_correlated(x[:, i])

          else:
              raise ValueError(f"Distribution '{namedist}' not supported.")

          # Update weights and fxixj for variable i
          w, fx_over_norm = update_weights(fx, hx, zf[:, i], zk_col)
          weight *= w
          fxixj *= fx_over_norm

      # Final correction with multivariate PDFs
      norm_multivarf = multivariate_normal(cov=matrix)
      phif = norm_multivarf.pdf(zf)
      norm_multivarh = multivariate_normal(cov=matrix)
      phih = norm_multivarh.pdf(zk)
      weight *= phif / phih
      fxixj *= phif

      return x, weight, fxixj 
  
  def var_rvs(self, ns, indexes_uncorrelated_xvar, nsigma=1.00, iprint=False):
    """

        Random variables generator for the Monte Carlo Simulation methods, only to uncorrelated variables

    """

    #Get only correlated variables
    xvar_uncorrelated = [self.reliability.xvar[i] for i in indexes_uncorrelated_xvar]
    nxvar_uncorrelated = len(xvar_uncorrelated)

    def fkapa(kapa, deltax, gsignal):
        fk = 1.00 + deltax ** 2 - gamma(1.00 + gsignal * 2.00 / kapa) / gamma(1.00 + gsignal * 1.00 / kapa) ** 2
        return fk
    
    def beta_limits(vars, mux, sigmax, q, r):
        a, b = vars
        eq1 = a + q / (q + r) * (b - a) - mux
        eq2 = ((q * r) / ((q + r) ** 2 * (q + r + 1))) ** (0.50) * (b - a) - sigmax
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

    #L = scipy.linalg.cholesky(self.reliability.correlation.Rz_rectify, lower=True)
    #Jzy = np.copy(L)

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
            print(self.reliability.xvar[i])
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
            
            ##Tem que arrumar no método de Bucher e Adaptive
            print('mufx', float(var['varmean']), self.reliability.xvarClassUncorrelated[i].mufx)
            print('sigmafx', float(var['varstd']), self.reliability.xvarClassUncorrelated[i].sigmafx)
            print('muhx', float(var['varhmean']), self.reliability.xvarClassUncorrelated[i].muhx)
            print('sigmahx', nsigma * sigmafx, self.reliability.xvarClassUncorrelated[i].sigmahx)
            
            x[:, i] = self.reliability.xvarClassUncorrelated[i].x_uncorrelated(ns)
            fx = self.reliability.xvarClassUncorrelated[i].fx(x[:, i])
            hx = self.reliability.xvarClassUncorrelated[i].hx(x[:, i])
            weight = weight * (fx / hx)
            fxixj = fxixj * fx 
        #
        # Uniform or constant distribution
        #
        
        elif namedist.lower() == 'uniform':    
            x[:, i] = self.reliability.xvarClassUncorrelated[i].x_uncorrelated(ns)
            fx = self.reliability.xvarClassUncorrelated[i].fx(x[:, i])
            hx = self.reliability.xvarClassUncorrelated[i].hx(x[:, i])
            weight = weight * (fx / hx)
            fxixj = fxixj * fx 
        #
        # Lognormal distribution
        #
        elif namedist.lower() == 'lognorm':
            
            x[:, i] = self.reliability.xvarClassUncorrelated[i].x_uncorrelated(ns)
            fx = self.reliability.xvarClassUncorrelated[i].fx_uncorrelated(x[:, i])
            hx = self.reliability.xvarClassUncorrelated[i].hx_uncorrelated(x[:, i])
            weight = weight * (fx / hx)
            fxixj = fxixj * fx 

        #
        # Gumbel distribution
        #
        elif namedist.lower() == 'gumbel':
            
            x[:, i] = self.reliability.xvarClassUncorrelated[i].x_uncorrelated(ns)
            fx = self.reliability.xvarClassUncorrelated[i].fx(x[:, i])
            hx = self.reliability.xvarClassUncorrelated[i].hx(x[:, i])
            weight = weight * (fx / hx)
            fxixj = fxixj * fx 

        #
        # Frechet distribution
        #
        elif namedist.lower() == 'frechet':
            
            # deltafx = sigmafx / mufx
            # kapa0 = 2.50
            # gsinal = -1.00
            # kapaf = scipy.optimize.newton(fkapa, kapa0, args=(deltafx, gsinal))
            # vfn = mufx / gamma(1.00 - 1.00 / kapaf)
            # deltahx = sigmahx / muhx
            # kapa0 = 2.50
            # gsinal = -1.00
            # kapah = scipy.optimize.newton(fkapa, kapa0, args=(deltahx, gsinal))
            # vhn = muhx / gamma(1.00 - 1.00 / kapah)
            
            
            x[:, i] = self.reliability.xvarClassUncorrelated[i].x_uncorrelated(ns)
            fx = self.reliability.xvarClassUncorrelated[i].fx_uncorrelated(x[:, i])
            hx = self.reliability.xvarClassUncorrelated[i].hx_uncorrelated(x[:, i])
            weight = weight * (fx / hx)
            fxixj = fxixj * fx 

        #
        #
        # Weibull distribution - minimum
        #
        elif namedist.lower() == 'weibull':
            # mufx = float(var['varmean'])
            # sigmafx = float(var['varstd'])
            # epsilon = float(var['varinf'])
            # muhx = float(var['varhmean'])
            # sigmahx = nsigma * sigmafx
            
            # deltafx = sigmafx / (mufx - epsilon)
            # kapa0 = 2.50
            # gsinal = 1.00
            # kapaf = scipy.optimize.newton(fkapa, kapa0, args=(deltafx, gsinal))
            # w1f = (mufx - epsilon) / gamma(1.00 + 1.00 / kapaf) + epsilon
            # deltahx = sigmahx / (muhx - epsilon)
            # kapa0 = 2.50
            # gsinal = 1.00
            
            # kapah = scipy.optimize.newton(fkapa, kapa0, args=(deltahx, gsinal))
            # w1h = (muhx - epsilon) / gamma(1.00 + 1.00 / kapah) + epsilon
            
            
            x[:, i] = self.reliability.xvarClassUncorrelated[i].x_uncorrelated(ns)
            fx = self.reliability.xvarClassUncorrelated[i].fx_uncorrelated(x[:, i])
            hx = self.reliability.xvarClassUncorrelated[i].hx_uncorrelated(x[:, i])
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
            muhx = float(var['varhmean'])
            sigmahx = nsigma * sigmafx
            loc = a
            scale = (b - a)
          
          ##Tem que arrumar no método de Bucher e Adaptive
            print('mufx', float(var['varmean']), self.reliability.xvarClassUncorrelated[i].mufx)
            print('sigmafx', float(var['varstd']), self.reliability.xvarClassUncorrelated[i].sigmafx)
            print('muhx', float(var['varhmean']), self.reliability.xvarClassUncorrelated[i].muhx)
            print('sigmahx', nsigma * sigmafx, self.reliability.xvarClassUncorrelated[i].sigmahx)
            
            
            
            
            # ah, bh =  fsolve(beta_limits, (1, 1), args= ( muhx, sigmahx, q, r))  
            # loch = ah
            # scaleh = (bh - ah)  
            
            
                  
            x[:, i] = self.reliability.xvarClassUncorrelated[i].x_uncorrelated(ns)
            fx = self.reliability.xvarClassUncorrelated[i].fx_uncorrelated(x[:, i])
            hx = self.reliability.xvarClassUncorrelated[i].hx_uncorrelated(x[:, i])
            weight = weight * (fx / hx)
            fxixj = fxixj * fx 

        #
        #
        # Gamma distribution
        #
        elif namedist.lower() == 'gamma':
            
            # mufx = float(var['varmean'])
            # sigmafx = float(var['varstd'])
            # deltafx = sigmafx / mufx
            # k = 1. / deltafx ** 2
            # v = k / mufx
            # a = k
            # loc = 0.00
            # scale = 1. / v
            # muhx = float(var['varhmean'])
            # sigmahx = nsigma * sigmafx
            # deltahx = sigmahx / muhx
            # kh = 1. / deltahx ** 2
            # vh = kh / muhx
            # ah = kh
            # loch = 0.00
            # scaleh = 1. / vh
            x[:, i] = self.reliability.xvarClassUncorrelated[i].x_uncorrelated(ns)
            fx = self.reliability.xvarClassUncorrelated[i].fx_uncorrelated(x[:, i])
            hx = self.reliability.xvarClassUncorrelated[i].hx_uncorrelated(x[:, i])
            weight = weight * (fx / hx)
            fxixj = fxixj * fx 
            

    return x, weight, fxixj
  

