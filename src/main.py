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
from scipy import optimize
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
import scipy.linalg
from scipy.special import gamma
import pandas as pd
import matplotlib.pyplot as plt
import time
from utils.distribution import createDistribution
from visualize import DataVisualize
from utils.validate.domain_types.validate_dvar import ValidateDvar
from utils.validate.domain_types.validate_gx import ValidateGx
from utils.validate.domain_types.validate_corrmatrix import ValidateCorrelationMatrix
from correlation import Correlation


class Reliability(Correlation):

    def __init__(self, xvar, dvar, gx, x0=None, Rz=None):
        
        #Simple setters
        self.dvar = dvar 
        self.fel = gx
        self.ndvar = len(dvar)
        self.nxvar = len(xvar)

        #Data processing - Complex setters
        self._xvar = self.set_xvar(xvar, x0)
        self._Rz = self.set_Rz(Rz, self._xvar, self.nxvar)
        self._d = self.set_d(self.dvar, self.ndvar)
        self._x0 = self.set_x0(x0, self.nxvar, self.xvar)

        super().__init__()        
           
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


    def set_Rz(self, Rz, xvar, nxvar):
      if Rz is None:
          return np.eye(self.nxvar)
      else:
        ValidateCorrelationMatrix(Rz)
        return self.nataf(Rz, xvar, nxvar)
        
    @nxvar.setter
    def nxvar(self, nxvar): 
      self._nxvar = nxvar
  
    @ndvar.setter
    def ndvar(self, ndvar): 
      self._ndvar = ndvar

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

    def form(self, iHLRF, toler=1.e-3, iprint=True):
        """

               Algorithm FORM-iHLRF. Normal equivalente transformation

        """
      #
        # FORM - First Order Reliability Method with improved HLRF (iHLRF)
        #
        #
        #
        # Penalty function m(y) for FORM-iHLRF algorithm
        #

        def mfunc(normy, g, c):
            my = 1. / 2. * normy ** 2 + c * np.abs(g)
            return my

        #
        #
        # Evaluation of parameter k for Frechet and Weibull distribution min
        #

        def fkapa(kapa, deltax, gsignal):
            fk = 1.00 + deltax ** 2 - gamma(1.00 + gsignal * 2.00 / kapa) / gamma(1.00 + gsignal * 1.00 / kapa) ** 2
            return fk

        #
        # Equivalent normal distribution parameters
        # xval = value of the variable x (scalar)
        # xpar1,xpar2,xpar3,xpar4 = parameters of the original pdf (scalars)
        # namedist = name of the x probability distribution ('string')
        #

        def normeqv(xval, xpar1, xpar2, xpar3, xpar4, namedist):

            #
            # Normal distribution
            #
            if namedist.lower() == 'gauss':
                mux = xpar1
                sigmax = xpar2
                muxneq = mux
                sigmaxneq = sigmax
            #
            # Uniform or constant distribution
            #
            elif namedist.lower() == 'uniform':
                a = xpar1
                b = xpar2
                c = (b - a)
                pdfx = 1. / c
                cdfx = (xval - a) / c
                zval = norm.ppf(cdfx)
                sigmaxneq = (norm.pdf(zval)) / pdfx
                muxneq = xval - zval * sigmaxneq
            #
            # Lognormal distribution
            #
            elif namedist.lower() == 'lognorm':
                mux = xpar1
                sigmax = xpar2
                zetax = np.sqrt(np.log(1. + (sigmax / mux) ** 2))
                lambdax = np.log(mux) - 0.50 * zetax ** 2
                sigmaxneq = zetax * xval
                muxneq = xval * (1. - np.log(xval) + lambdax)
            #
            # Gumbel distribution
            #
            elif namedist.lower() == 'gumbel':
                mux = xpar1
                sigmax = xpar2
                alphan = (np.pi / np.sqrt(6.00)) / (sigmax)
                un = mux - np.euler_gamma / alphan
                cdfx = np.exp(-np.exp(-alphan * (xval - un)))
                pdfx = alphan * np.exp(-alphan * (xval - un)) * cdfx
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq
            #
            #
            # Frechet distribution
            #
            elif namedist.lower() == 'frechet':
                mux = xpar1
                sigmax = xpar2
                deltax = sigmax / mux
                kapa0 = 2.50
                gsignal = -1.00
                kapa = scipy.optimize.newton(fkapa, kapa0, args=(deltax, gsignal))
                vn = mux / gamma(1.00 - 1.00 / kapa)
                cdfx = np.exp(-(vn / xval) ** kapa)
                pdfx = kapa / vn * (vn / xval) ** (kapa + 1) * np.exp(-(vn / xval) ** kapa)
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq
            #
            #
            # Weibull distribution minimum
            #
            elif namedist.lower() == 'weibull':
                mux = xpar1
                sigmax = xpar2
                epsilon = xpar3
                deltax = sigmax / (mux - epsilon)
                kapa0 = 2.50
                gsignal = 1.00
                kapa = scipy.optimize.newton(fkapa, kapa0, args=(deltax, gsignal))
                w1 = (mux - epsilon) / gamma(1.00 + 1.00 / kapa) + epsilon
                y1 = (xval - epsilon) / (w1 - epsilon)
                pdfx = weibull_min.pdf(y1, kapa) / (w1 - epsilon)
                cdfx = weibull_min.cdf(y1, kapa)
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq
            #
            #
            # Beta distribution
            #
            elif namedist.lower() == 'beta':
                a = xpar1
                b = xpar2
                q = xpar3
                r = xpar4
                loc = a
                scale = (b - a)
                pdfx = beta_dist.pdf(xval, q, r, loc, scale)
                cdfx = beta_dist.cdf(xval, q, r, loc, scale)
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq


            #
            #
            # Gamma distribution
            #
            elif namedist.lower() == 'gamma':
                mux = xpar1
                sigmax = xpar2
                delta = sigmax / mux
                k = 1. / delta ** 2
                v = k / mux
                a = k
                loc = 0.00
                scale = 1. / v
                pdfx = gamma_dist.pdf(xval, a, loc, scale)
                cdfx = gamma_dist.cdf(xval, a, loc, scale)
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq

            return muxneq, sigmaxneq
        
        #
        #
        # Data input
        #
        # Number of variables of the problem

        # Equivalent normal mean and standard deviation of the variables
        muxneqk = np.zeros(self.nxvar)
        sigmaxneqk = np.zeros(self.nxvar)
        namevar = []
        dist = []
        mux0 = []
        sigmax0 = []
        par1 = []
        par2 = []
        par3 = []
        par4 = []
        #
        # Original mean and standard deviation of the variables x
        #

        i = -1
        for var in self.xvar:
            i += 1
            # Names of the random variables x
            namevar.append(str(var['varname']))
            # Names of the probability density functions of the variables x
            dist.append(str(var['vardist']))
            # Mean value of the random variables x
            mux0.append(float(var['varmean']))
            # Standard deviation of the random variables x
            if var['varstd'] == 0.00:
                sigmax0.append(float(var['varcov']) * float(var['varmean']))
            else:
                sigmax0.append(float(var['varstd']))
            # Parameter1
            if 'parameter1' in var:
                par1.append(float(var['parameter1']))
            else:
                par1.append(0.00)
            # Parameter2
            if 'parameter2' in var:
                par2.append(float(var['parameter2']))
            else:
                par2.append(0.00)
            # Parameter3
            if 'parameter3' in var:
                par3.append(float(var['parameter3']))
            else:
                par3.append(0.00)
            # Parameter4
            if 'parameter4' in var:
                par4.append(float(var['parameter4']))
            else:
                par4.append(0.00)
            
           
        #
        # Conversion to array format
        #
        mux0 = np.array(mux0)
        sigmax0 = np.array(sigmax0)
        par1 = np.array(par1)
        par2 = np.array(par2)
        par3 = np.array(par3)
        par4 = np.array(par4)
        #
        #   Algorithm FORM-HLRF: Beck, 2019, pag. 101.
        #
        #
        # Step 1 - Determination of equivalent correlation coefficients and
        #          Jacobian matrices Jxz and Jzx
        #
        Imatrix = np.eye(self.nxvar)
        #
        # Correlation matrix is self.corrmatrix
        #
        if iprint:
            print('Correlation Matrix after Nataf correction:')
            print(self.Rz)
        #
        # Cholesky decomposition of the correlation matrix
        #
        L = scipy.linalg.cholesky(self.Rz, lower=True)
        Jzy = np.copy(L)
        Jyz = np.linalg.inv(L)
        #
        # Step 2 - Initialize de xk value with mux0
        #
        # Initialization of the variable yk1
        # Jacobian matrices of x==>z and z==>y transformations
        D = sigmax0 * Imatrix
        Jzx = np.linalg.inv(D)
        Jyx = np.dot(Jyz, Jzx)
        Jxz = np.copy(D)
        Jxy = np.dot(Jxz, Jzy)
        yk1 = np.zeros(self.nxvar)
    #    xk1 = mux0 + Jxy.dot(yk1)
        xk1 = np.copy(self.x0)
        #
        # Error tolerance for yk and g(x)
        epsilon = toler
        delta = toler * np.abs(self.fel(xk1, self.d))
        # Initial values for errors and iteration counters
        erro1 = 1000.00
        erro2 = 1000.00
        kiter = 0
        # Value of dx increment for the evaluation of the derivatives
        eps = 1.e-6
        #
        while (erro1 > epsilon or erro2 > delta) and kiter < 100:
            #
            kiter += 1
            xk = np.copy(xk1)
            #
            # Calculation of the equivalent normal distribution parameters for xk
            #
            for i in range(self.nxvar):
                xval = xk[i]
                mux = mux0[i]
                sigmax = sigmax0[i]
                namedist = dist[i]
                xpar1 = mux
                xpar2 = sigmax
                xpar3 = par3[i]
                xpar4 = par4[i]
                if dist[i] == 'beta':
                    xpar1 = par1[i]
                    xpar2 = par2[i]
                if dist[i] == 'uniform':
                    xpar1 = par1[i]
                    xpar2 = par2[i]

                muxneqk[i], sigmaxneqk[i] = normeqv(xval, xpar1, xpar2, xpar3, xpar4, namedist)
            #
            # Step 3 - Update of the Jacobian matrices Jyx and Jxy
            #
            Dneq = sigmaxneqk * Imatrix
            Jzx = np.linalg.inv(Dneq)
            Jyx = np.dot(Jyz, Jzx)
            Jxz = np.copy(Dneq)
            Jxy = np.dot(Jxz, Jzy)
            #
            #  Step 4 - Transformation from xk to yk
            #
            yk = Jyx.dot(xk - muxneqk)
            normyk = np.linalg.norm(yk)
            beta = np.linalg.norm(yk)

            #
            #  Step 5 - Evaluation of g(xk)
            #
            gxk = self.fel(xk, self.d)

            #
            # Step 6 - Evaluation of the gradients of g(x) in relation to yk
            #
            #
            # a. Calculation of the partial derivatives of g(x) in relation to xk
            #
            gradxk = optimize.approx_fprime(xk, self.fel, eps, self.d)
            #
            # b. Calculation of the partial derivatives of g(x) in relation to yk
            #
            gradyk = np.transpose(Jxy).dot(gradxk)
            normgradyk = np.linalg.norm(gradyk)
            #
            # c. Calculation of the direction cosines for xk
            #
            # Direction cosines
            alpha = gradyk / normgradyk

            #
            # Step 7. Vector yk updating to yk+1 by HLRF algorithm
            #
            dk = ((np.dot(gradyk, yk) - gxk) / normgradyk ** 2) * gradyk - yk
            lambdak = 1.00
            yk1 = yk + lambdak * dk
            #
            # Parameters of iHLRF method
            #
            if iHLRF:
                gamma0 = 2.0
                a = 0.1
                # a = 0.5
                b = 0.5
                #
                gyk = gxk
                normyk = np.linalg.norm(yk)
                normyk1 = np.linalg.norm(yk1)
                c1 = normyk / normgradyk
                #
                if erro2 > delta:
                    c2 = 0.5 * normyk1 ** 2 / np.abs(gyk)
                    ck = gamma0 * np.max([c1, c2])
                else:
                    ck = gamma0 * c1
                #
                k = -1
                f1 = 1.00
                f2 = 0.00
                while f1 > f2 and k < 10:
                    k += 1
                    lambdak = b ** k
                    yk1 = yk + lambdak * dk
                    xk1 = muxneqk + Jxy.dot(yk1)
                    gyk1 = self.fel(xk1, self.d)
                    normyk1 = np.linalg.norm(yk1)
                    f1 = mfunc(normyk1, gyk1, ck) - mfunc(normyk, gyk, ck)
                    gradm = yk + ck * gradyk * np.sign(gyk)
                    normgradm = np.linalg.norm(gradm)
                    f2 = a * lambdak * np.dot(gradm, dk)
            #        f2=-a*lambdak*normgradm**2 # Beck pg. 85: It does not work!!
            #        res=np.array([k,ck,lambdak,gxk,gyk1,f1,f2])
            #        print(res)
            #
            yk1 = yk + lambdak * dk

            #
            # Step 8. Transformation from yk+1 to xk+1
            #
            xk1 = muxneqk + Jxy.dot(yk1)

            #
            # Step 9. Convergence test for yk and g(x)
            #
            prod = normgradyk * normyk
            # Evaluation of the error in the yk1 vector
            if np.abs(prod) > eps:
                erro1 = 1. - np.abs(np.dot(gradyk, yk) / (normgradyk * normyk))
            else:
                erro1 = 1000.00
            # Evaluation of the error in the limit state function g(x)
            erro2 = np.abs(gxk)
            # Printing of the updated values
            if iprint:
                print('\nIteration number = {0:d} g(x) ={1:0.5e} erro1 ={2:0.5e} Beta ={3:0.4f}'
                      .format(kiter, gxk, erro1, beta))
                datadict = {'xvar': namevar, 'prob_dist': dist, 'mux': muxneqk, 'sigmax': sigmaxneqk,
                            'xk': xk, 'yk': yk, 'alpha': alpha}
                data = pd.DataFrame(datadict)
                print(data)
        #
        pf = norm.cdf(-beta)
        if iprint:
            print('\nProbability of Failure Pf = {0:0.4e}'.format(pf))

        return {
          "beta": beta,
          "yk":yk,
          "xk": xk,
          "alpha": alpha,
          "normgradyk": normgradyk,
          "sigmaxneqk": sigmaxneqk
          } 

    def form2(self, iHLRF, toler=1.e-3, iprint=True):
        """

            Algorithm FORM-iHLRF. Direct mapping to standard Gaussian space

        """
      #
        # FORM - First Order Reliability Method with improved HLRF (iHLRF)
        #
        #
        #
        # Penalty function m(y) for FORM-iHLRF algorithm
        #

        def mfunc(normy, g, c):
            my = 1. / 2. * normy ** 2 + c * np.abs(g)
            return my

        #
        #
        # Evaluation of parameter k for Frechet and Weibull distribution minimum
        #

        def fkapa(kapa, deltax, gsignal):
            fk = 1.00 + deltax ** 2 - gamma(1.00 + gsignal * 2.00 / kapa) / gamma(1.00 + gsignal * 1.00 / kapa) ** 2
            return fk

        #
        # Equivalent normal distribution parameters
        # xval = value of the variable x (scalar)
        # xpar1,xpar2,xpar3,xpar4 = parameters of the original pdf (scalars)
        # namedist = name of the x probability distribution ('string')
        # zval = equivalente normal variabel correlated
        #
        def normeqv(xval, xpar1, xpar2, xpar3, xpar4, namedist):

            #
            # Normal distribution
            #
            if namedist.lower() == 'gauss':
                mux = xpar1
                sigmax = xpar2
                muxneq = mux
                sigmaxneq = sigmax
            #
            # Uniform or constant distribution
            #
            elif namedist.lower() == 'uniform':
                epsilon = 1.e-8
                a = xpar1
                b = xpar2
                c = (b - a)
                if xval<=a: xval = a + epsilon
                if xval>=b: xval = b - epsilon
                pdfx = 1. / c
                cdfx = (xval - a) / c
                zval = norm.ppf(cdfx)
                sigmaxneq = (norm.pdf(zval)) / pdfx
                muxneq = xval - zval * sigmaxneq
            #
            # Lognormal distribution
            #
            elif namedist.lower() == 'lognorm':
                mux = xpar1
                sigmax = xpar2
                zetax = np.sqrt(np.log(1. + (sigmax / mux) ** 2))
                lambdax = np.log(mux) - 0.50 * zetax ** 2
                sigmaxneq = zetax * xval
                muxneq = xval * (1. - np.log(xval) + lambdax)
            #
            # Gumbel distribution
            #
            elif namedist.lower() == 'gumbel':
                mux = xpar1
                sigmax = xpar2
                alphan = (np.pi / np.sqrt(6.00)) / (sigmax)
                un = mux - np.euler_gamma / alphan
                cdfx = np.exp(-np.exp(-alphan * (xval - un)))
                pdfx = alphan * np.exp(-alphan * (xval - un)) * cdfx
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq
            #
            #
            # Frechet distribution
            #
            elif namedist.lower() == 'frechet':
                mux = xpar1
                sigmax = xpar2
                deltax = sigmax / mux
                kapa0 = 2.50
                gsignal = -1.00
                kapa = scipy.optimize.newton(fkapa, kapa0, args=(deltax, gsignal))
                vn = mux / gamma(1.00 - 1.00 / kapa)
                cdfx = np.exp(-(vn / xval) ** kapa)
                pdfx = kapa / vn * (vn / xval) ** (kapa + 1) * np.exp(-(vn / xval) ** kapa)
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq
            #
            #
            # Weibull distribution minimum
            #
            elif namedist.lower() == 'weibull':
                mux = xpar1
                sigmax = xpar2
                epsilon = xpar3
                deltax = sigmax / (mux - epsilon)
                kapa0 = 2.50
                gsignal = 1.00
                kapa = scipy.optimize.newton(fkapa, kapa0, args=(deltax, gsignal))
                w1 = (mux - epsilon) / gamma(1.00 + 1.00 / kapa) + epsilon
                y1 = (xval - epsilon) / (w1 - epsilon)
                pdfx = weibull_min.pdf(y1, kapa) / (w1 - epsilon)
                cdfx = weibull_min.cdf(y1, kapa)
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq
            #
            #
            # Beta distribution
            #
            elif namedist.lower() == 'beta':
                epsilon = 1.e-8
                a = xpar1
                b = xpar2
                q = xpar3
                r = xpar4
                loc = a
                scale = (b - a)
                if xval<=a: xval = a + epsilon
                if xval>=b: xval = b - epsilon
                pdfx = beta_dist.pdf(xval, q, r, loc, scale)
                cdfx = beta_dist.cdf(xval, q, r, loc, scale)
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq

            #
            #
            # Gamma distribution
            #
            elif namedist.lower() == 'gamma':
                mux = xpar1
                sigmax = xpar2
                deltax = sigmax / mux
                k = 1. / deltax
                v = k / mux
                a = k
                loc = 0.00
                scale = 1. / v
                pdfx = gamma_dist.pdf(xval, a, loc, scale)
                cdfx = gamma_dist.cdf(xval, a, loc, scale)
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq

            return muxneq, sigmaxneq
        
        def direct_mapping(xval, xpar1, xpar2, xpar3, xpar4, namedist):

            #
            # Normal distribution. Direct mapping to standard Gaussian space
            #
            if namedist.lower() == 'gauss':
                mux = xpar1
                sigmax = xpar2
                cdfx = norm.cdf(xval,mux,sigmax)
                
            #
            # Uniform or constant distribution. Direct mapping to standard Gaussian space
            #
            elif namedist.lower() == 'uniform':
                epsilon = 1.e-8
                a = xpar1
                b = xpar2
                c = (b - a)
                if xval<=a: xval = a + epsilon
                if xval>=b: xval = b - epsilon
                cdfx = uniform.cdf(xval,a,c)
                
                
            #
            # Lognormal distribution. Direct mapping to standard Gaussian space
            #
            elif namedist.lower() == 'lognorm':
                mux = xpar1
                sigmax = xpar2
                zetax = np.sqrt(np.log(1. + (sigmax / mux) ** 2))
                lambdax = np.log(mux) - 0.50 * zetax ** 2
                cdfx = lognorm.cdf(xval/np.exp(lambdax), zetax)
                
            #
            # Gumbel distribution. Direct mapping to standard Gaussian space
            #
            elif namedist.lower() == 'gumbel':
                mux = xpar1
                sigmax = xpar2
                alphan = (np.pi / np.sqrt(6.00)) / (sigmax)
                un = mux - np.euler_gamma / alphan
                betan = 1.00 / alphan
                cdfx = gumbel_r.cdf(xval, un, betan)
                
            #
            #
            # Frechet distribution. Direct mapping to standard Gaussian space
            #
            elif namedist.lower() == 'frechet':
                mux = xpar1
                sigmax = xpar2
                deltax = sigmax / mux
                kapa0 = 2.50
                gsignal = -1.00
                kapa = scipy.optimize.newton(fkapa, kapa0, args=(deltax, gsignal))
                vn = mux / gamma(1.00 - 1.00 / kapa)
                yn = xval / vn
                cdfx = invweibull.cdf(yn, kapa)
                
            #
            #
            # Weibull distribution minimum. Direct mapping to standard Gaussian space
            #
            elif namedist.lower() == 'weibull':
                mux = xpar1
                sigmax = xpar2
                epsilon = xpar3
                deltax = sigmax / (mux - epsilon)
                kapa0 = 2.50
                gsignal = 1.00
                kapa = scipy.optimize.newton(fkapa, kapa0, args=(deltax, gsignal))
                w1 = (mux - epsilon) / gamma(1.00 + 1.00 / kapa) + epsilon
                yn = (xval - epsilon) / (w1 - epsilon)
                cdfx = weibull_min.cdf(yn, kapa)
                
            #
            #
            # Beta distribution. Direct mapping to standard Gaussian space
            #
            elif namedist.lower() == 'beta':
                epsilon = 1.e-8
                a = xpar1
                b = xpar2
                q = xpar3
                r = xpar4
                loc = a
                scale = (b - a)
                if xval<=a: xval = a + epsilon
                if xval>=b: xval = b - epsilon
                cdfx = beta_dist.cdf(xval, q, r, loc, scale)

            #
            #
            # Gamma distribution. Direct mapping to standard Gaussian space
            #
            elif namedist.lower() == 'gamma':
                mux = xpar1
                sigmax = xpar2
                delta = sigmax / mux
                k = 1. / delta ** 2
                v = k / mux
                a = k
                loc = 0.00
                scale = 1. / v
                cdfx = gamma_dist.cdf(xval, a, loc, scale)
                

            return cdfx
        
        
        def inverse_mapping(zval, xpar1, xpar2, xpar3, xpar4, namedist):

            #
            # Normal distribution. Inverse mapping from standard Gaussian space
            #
            if namedist.lower() == 'gauss':
                mux = xpar1
                sigmax = xpar2
                xval = norm.ppf(cdfx,mux,sigmax)
            #
            # Uniform or constant distribution. Inverse mapping from standard Gaussian space
            #
            elif namedist.lower() == 'uniform':
                a = xpar1
                b = xpar2
                c = (b - a)
                xval = uniform.ppf(cdfx,a,c)
            #
            # Lognormal distribution. Inverse mapping from standard Gaussian space
            #
            elif namedist.lower() == 'lognorm':
                mux = xpar1
                sigmax = xpar2
                zetax = np.sqrt(np.log(1. + (sigmax / mux) ** 2))
                lambdax = np.log(mux) - 0.50 * zetax ** 2
                xval = lognorm.ppf(cdfx, zetax)*np.exp(lambdax)
            #
            # Gumbel distribution. Inverse mapping from standard Gaussian space
            #
            elif namedist.lower() == 'gumbel':
                mux = xpar1
                sigmax = xpar2
                alphan = (np.pi / np.sqrt(6.00)) / (sigmax)
                un = mux - np.euler_gamma / alphan
                xval = un - 1./alphan * np.log(np.log(1. / cdfx))
            #
            #
            # Frechet distribution. Inverse mapping from standard Gaussian space
            #
            elif namedist.lower() == 'frechet':
                mux = xpar1
                sigmax = xpar2
                deltax = sigmax / mux
                kapa0 = 2.50
                gsignal = -1.00
                kapa = scipy.optimize.newton(fkapa, kapa0, args=(deltax, gsignal))
                vn = mux / gamma(1.00 - 1.00 / kapa)
                xval = vn / (np.log(1. / cdfx)) ** (1. / kapa)
            #
            #
            # Weibull distribution minimum. Inverse mapping from standard Gaussian space
            #
            elif namedist.lower() == 'weibull':
                mux = xpar1
                sigmax = xpar2
                epsilon = xpar3
                deltax = sigmax / (mux - epsilon)
                kapa0 = 2.50
                gsignal = 1.00
                kapa = scipy.optimize.newton(fkapa, kapa0, args=(deltax, gsignal))
                w1 = (mux - epsilon) / gamma(1.00 + 1.00 / kapa) + epsilon
                xval = (w1 - epsilon) * (np.log(1./(1. - cdfx))) ** (1. / kapa) + epsilon
            #
            #
            # Beta distribution. Inverse mapping from standard Gaussian space
            #
            elif namedist.lower() == 'beta':
                a = xpar1
                b = xpar2
                q = xpar3
                r = xpar4
                loc = a
                scale = (b - a)
                xval = beta_dist.ppf(cdfx, q, r, loc, scale)

            #
            #
            # Gamma distribution. Inverse mapping from standard Gaussian space
            #
            elif namedist.lower() == 'gamma':
                mux = xpar1
                sigmax = xpar2
                delta = sigmax / mux
                k = 1. / delta ** 2
                v = k / mux
                a = k
                loc = 0.00
                scale = 1. / v
                xval = gamma_dist.ppf(cdfx, a, loc, scale)
                

            return xval
        #
        #
        # Data input
        #
        # Number of variables of the problem

        # Equivalent normal mean and standard deviation of the variables
        muxneqk = np.zeros(self.nxvar)
        sigmaxneqk = np.zeros(self.nxvar)
        namevar = []
        dist = []
        mux0 = []
        sigmax0 = []
        par1 = []
        par2 = []
        par3 = []
        par4 = []
        #
        # Original mean and standard deviation of the variables x
        #

        i = -1
        for var in self.xvar:
            i += 1
            # Names of the random variables x
            namevar.append(str(var['varname']))
            # Names of the probability density functions of the variables x
            dist.append(str(var['vardist']))
            # Mean value of the random variables x
            mux0.append(float(var['varmean']))
            # Standard deviation of the random variables x
            if var['varstd'] == 0.00:
                sigmax0.append(float(var['varcov']) * float(var['varmean']))
            else:
                sigmax0.append(float(var['varstd']))
            # Parameter1
            if 'parameter1' in var:
                par1.append(float(var['parameter1']))
            else:
                par1.append(0.00)
            # Parameter2
            if 'parameter2' in var:
                par2.append(float(var['parameter2']))
            else:
                par2.append(0.00)
            # Parameter3
            if 'parameter3' in var:
                par3.append(float(var['parameter3']))
            else:
                par3.append(0.00)
            # Parameter4
            if 'parameter4' in var:
                par4.append(float(var['parameter4']))
            else:
                par4.append(0.00)
            
           
        #
        # Conversion to array format
        #
        mux0 = np.array(mux0)
        sigmax0 = np.array(sigmax0)
        par1 = np.array(par1)
        par2 = np.array(par2)
        par3 = np.array(par3)
        par4 = np.array(par4)
        #
        #   Algorithm FORM-HLRF: Beck, 2019, pag. 101.
        #
        #
        # Step 1 - Determination of equivalent correlation coefficients and
        #          Jacobian matrices Jxz and Jzx
        #
        Imatrix = np.eye(self.nxvar)
        #
        # Correlation matrix is self.corrmatrix
        #
        if iprint:
            print('Correlation Matrix after Nataf correction:')
            print(self.Rz)
        #
        # Cholesky decomposition of the correlation matrix
        #
        L = scipy.linalg.cholesky(self.Rz, lower=True)
        Jzy = np.copy(L)
        Jyz = np.linalg.inv(L)
        #
        # Step 2 - Initialize de xk value with mux0
        #
        # Initialization of the variable yk1
        # Jacobian matrices of x==>z and z==>y transformations
        D = sigmax0 * Imatrix
        Jzx = np.linalg.inv(D)
        Jyx = np.dot(Jyz, Jzx)
        Jxz = np.copy(D)
        Jxy = np.dot(Jxz, Jzy)
        yk1 = np.zeros(self.nxvar)
        zk = np.zeros(self.nxvar)
        zk1 = np.zeros(self.nxvar)
        uk = np.zeros(self.nxvar)
        uk1 = np.zeros(self.nxvar)
        xk1 = np.copy(self.x0)
        #
        # Error tolerance for yk and g(x)
        epsilon = toler
        delta = toler * np.abs(self.fel(xk1, self.d))
        # Initial values for errors and iteration counters
        erro1 = 1000.00
        erro2 = 1000.00
        kiter = 0
        # Value of dx increment for the evaluation of the derivatives
        eps = 1.e-6
        #
        while (erro1 > epsilon or erro2 > delta) and kiter < 100:
            #
            kiter += 1
            xk = np.copy(xk1)
            #
            # Step 2 - Calculation of equivalent normal variables by direct mapping to standard normal space (correlated)
            #
            for i in range(self.nxvar):
                xval = xk[i]
                mux = mux0[i]
                sigmax = sigmax0[i]
                namedist = dist[i]
                xpar1 = mux
                xpar2 = sigmax
                xpar3 = par3[i]
                xpar4 = par4[i]
                if dist[i] == 'beta':
                    xpar1 = par1[i]
                    xpar2 = par2[i]
                if dist[i] == 'uniform':
                    xpar1 = par1[i]
                    xpar2 = par2[i]

                uk[i] = direct_mapping(xval, xpar1, xpar2, xpar3, xpar4, namedist)

                muxneqk[i], sigmaxneqk[i] = normeqv(xval, xpar1, xpar2, xpar3, xpar4, namedist)
            #
            # Step 3 - Update of the Jacobian matrices Jyx and Jxy
            #
            Dneq = sigmaxneqk * Imatrix
            Jzx = np.linalg.inv(Dneq)
            Jyx = np.dot(Jyz, Jzx)
            Jxz = np.copy(Dneq)
            Jxy = np.dot(Jxz, Jzy)
            
            
            
            
            #
            #  Step 4 - Transformation from zk to yk
            #
            zk = norm.ppf(uk)
            yk = Jyz.dot(zk)
            normyk = np.linalg.norm(yk)
            beta = np.linalg.norm(yk)

            #
            #  Step 5 - Evaluation of g(xk)
            #
            gxk = self.fel(xk, self.d)

            #
            # Step 6 - Evaluation of the gradients of g(x) in relation to yk
            #
            #
            # a. Calculation of the partial derivatives of g(x) in relation to xk
            #
            gradxk = optimize.approx_fprime(xk, self.fel, eps, self.d)
            #
            # b. Calculation of the partial derivatives of g(x) in relation to yk
            #
            gradyk = np.transpose(Jxy).dot(gradxk)
            normgradyk = np.linalg.norm(gradyk)
            #
            # c. Calculation of the direction cosines for xk
            #
            # Direction cosines
            alpha = gradyk / normgradyk

            #
            # Step 7. Vector yk updating to yk+1 by HLRF algorithm
            #
            dk = ((np.dot(gradyk, yk) - gxk) / normgradyk ** 2) * gradyk - yk
            lambdak = 1.00
            yk1 = yk + lambdak * dk
            #
            # Parameters of iHLRF method
            #
            if iHLRF:
                gamma0 = 2.0
                a = 0.1
                # a = 0.5
                b = 0.5
                #
                gyk = gxk
                normyk = np.linalg.norm(yk)
                normyk1 = np.linalg.norm(yk1)
                c1 = normyk / normgradyk
                #
                if erro2 > delta:
                    c2 = 0.5 * normyk1 ** 2 / np.abs(gyk)
                    ck = gamma0 * np.max([c1, c2])
                else:
                    ck = gamma0 * c1
                #
                k = -1
                f1 = 1.00
                f2 = 0.00
                while f1 > f2 and k < 10:
                    k += 1
                    lambdak = b ** k
                    yk1 = yk + lambdak * dk
                    #
                    #  Transformation from yk+1 to xk+1
                    #
                    zk1 = Jzy.dot(yk1)
                    uk1 = norm.cdf(zk1)
                    for i in range(self.nxvar):
                        cdfx = uk1[i]
                        mux = mux0[i]
                        sigmax = sigmax0[i]
                        namedist = dist[i]
                        xpar1 = mux
                        xpar2 = sigmax
                        xpar3 = par3[i]
                        xpar4 = par4[i]
                        if dist[i] == 'beta':
                            xpar1 = par1[i]
                            xpar2 = par2[i]
                        if dist[i] == 'uniform':
                            xpar1 = par1[i]
                            xpar2 = par2[i]

                        xk1[i] = inverse_mapping(cdfx, xpar1, xpar2, xpar3, xpar4, namedist)

                    gyk1 = self.fel(xk1, self.d)
                    normyk1 = np.linalg.norm(yk1)
                    f1 = mfunc(normyk1, gyk1, ck) - mfunc(normyk, gyk, ck)
                    gradm = yk + ck * gradyk * np.sign(gyk)
                    normgradm = np.linalg.norm(gradm)
                    f2 = a * lambdak * np.dot(gradm, dk)
            
            #
            yk1 = yk + lambdak * dk

            #
            # Step 8. Transformation from yk+1 to xk+1
            #
            zk1 = Jzy.dot(yk1)
            uk1 = norm.cdf(zk1)
            for i in range(self.nxvar):
                cdfx = uk1[i]
                mux = mux0[i]
                sigmax = sigmax0[i]
                namedist = dist[i]
                xpar1 = mux
                xpar2 = sigmax
                xpar3 = par3[i]
                xpar4 = par4[i]
                if dist[i] == 'beta':
                    xpar1 = par1[i]
                    xpar2 = par2[i]
                if dist[i] == 'uniform':
                    xpar1 = par1[i]
                    xpar2 = par2[i]

                xk1[i] = inverse_mapping(cdfx, xpar1, xpar2, xpar3, xpar4, namedist)

            #
            # Step 9. Convergence test for yk and g(x)
            #
            prod = normgradyk * normyk
            # Evaluation of the error in the yk1 vector
            if np.abs(prod) > eps:
                erro1 = 1. - np.abs(np.dot(gradyk, yk) / (normgradyk * normyk))
            else:
                erro1 = 1000.00
            # Evaluation of the error in the limit state function g(x)
            erro2 = np.abs(gxk)
            # Printing of the updated values
            if iprint:
                print('\nIteration number = {0:d} g(x) ={1:0.5e} erro1 ={2:0.5e} Beta ={3:0.4f}'
                      .format(kiter, gxk, erro1, beta))
                datadict = {'xvar': namevar, 'prob_dist': dist, 'mux': muxneqk, 'sigmax': sigmaxneqk,
                            'xk': xk, 'yk': yk, 'alpha': alpha}
                data = pd.DataFrame(datadict)
                print(data)
        #
        pf = norm.cdf(-beta)
        if iprint:
            print('\nProbability of Failure Pf = {0:0.4e}'.format(pf))
   
        return {
            "beta": beta,
            "xk": xk,
            "alpha": alpha,
            "normgradyk": normgradyk,
            } 

    def sorm(self, iHLRF=True, toler=1.e-6, iprint=True):
        """
        Second order reliability method = SORM

        """

        #
        # GRAM-SCHMIDT transformation
        #
        def gramschmidt(A, nxvar):
            rk = np.zeros(nxvar)
            rj = np.zeros(nxvar)
            rk0 = np.zeros(nxvar)
            #
            R = np.zeros((nxvar, nxvar))
            R[nxvar - 1, :] = A[nxvar - 1, :].copy()
            for k in range(nxvar - 2, -1, -1):
                rk0 = A[k, :].copy()
                rk0projection = np.zeros(nxvar)
                for j in range(nxvar - 1, k, -1):
                    rj = R[j, :].copy()
                    projection = (rj.dot(rk0)) / (rj.dot(rj))
                    rk0projection = rk0projection + projection * rj
                rk = rk0 - rk0projection
                R[k, :] = rk.copy()
            for i in range(nxvar):
                R[i, :] = R[i, :] / np.linalg.norm(R[i, :])
            #
            return R

        #
        #
        # Function to calculate the second order derivative: d2g/dxidxj
        #
        def second_order_derivative(x, i, j):
            epsilon = 1.e-4  # tolerance for the increments
            h1 = epsilon  # increments: h1 and h2, when i is not equal to j
            h2 = epsilon  # different increments can be adopted
            h = epsilon  # increment h
            a = x[i]  # reference value for x[i]
            b = x[j]  # reference value for x[j]
            #
            # Code: gmn where m and n are equal to:
            # Index 0 = no increment is applied to the variables i and j
            # Index 1 = a decrement equal to -h is applied to the variable i (or j)
            # Index 2 = an incremente equal to +h is applied to the variable i (or j)
            #
            if i == j:
                x0 = np.copy(x)
                x0[i] = a - h
                g10 = self.fel(x0, self.d)
                x0[i] = a
                g00 = self.fel(x0, self.d)
                x0[i] = a + h
                g20 = self.fel(x0, self.d)
                d2g = (g10 - 2. * g00 + g20) / h ** 2  # second order derivative: d2g/dxi2
            else:
                x0 = np.copy(x)
                x0[i] = a + h1
                x0[j] = b + h2
                g22 = self.fel(x0, self.d)
                x0[i] = a + h1
                x0[j] = b - h2
                g21 = self.fel(x0, self.d)
                x0[i] = a - h1
                x0[j] = b + h2
                g12 = self.fel(x0, self.d)
                x0[i] = a - h1
                x0[j] = b - h2
                g11 = self.fel(x0, self.d)
                d2g = (g22 - g21 - g12 + g11) / (4. * h1 * h2)  # second order derivative: d2g/dxidxj
            #
            return d2g

        #
        # First run FORM-iHLRF algorithm
        #
        n = self.nxvar
        xk = np.zeros(n)
        yk = np.zeros(n)
        gradxk = np.zeros(n)
        alpha = np.zeros(n)
        beta = 0.00
        kiter = 0
        erro1 = 0.00

        beta, xk, alpha, normgradyk, sigmaxneqk = self.form(iHLRF, toler, iprint)
        #
        # Formulation of Second Order Reliability Method - SORM
        #
        #
        # Failure probability calculation
        #
        pfform = norm.cdf(-beta)
        #
        # Calculation of the Hessian Matrix
        #
        bmatrix = np.zeros((n, n))
        dmatrix = np.zeros((n, n))
        amatrix = np.eye(n)
        hmatrix = np.zeros((n, n))

        np.set_printoptions(precision=4)

        #
        # Calculation of the Hessian matrix D: d2g/dyidyj
        #
        for i in range(n):
            for j in range(n):
                dmatrix[i, j] = second_order_derivative(xk, i, j) * sigmaxneqk[i] * sigmaxneqk[j]


        #
        # Calculation of the matrix B
        #
        bmatrix = 1. / normgradyk * dmatrix

        #
        # Calculation of the orthogonal matrix H
        #
        amatrix[n - 1, :] = alpha.copy()
        #
        hmatrix = gramschmidt(amatrix, n)


        #
        # Calculation of the curvature matrix K
        #
        kmatrix = hmatrix.dot(bmatrix.dot(hmatrix.T))

        #
        # Calculation of the failure probability using SORM Breitung equation
        #
        factor = 1.00
        for i in range(n - 1):
            factor = factor * 1. / np.sqrt(1.00 + beta * kmatrix[i, i])
        pfsorm = pfform * factor
        betasorm = -norm.ppf(pfsorm)
        #
        # Print the result
        #
        if iprint:
            print('\nSORM results:')
            print('\nHessian matrix:')
            print(dmatrix)
            print('\nNorm of the gradient of g(y) =', normgradyk)
            print('\nB matrix:')
            print(bmatrix)
            print('\nH matrix:')
            print(hmatrix)
            print('\nK = curvatures matrix:')
            print(kmatrix)
            print('\npfFORM =', pfform)
            print('\nBetaFORM =', beta)
            print('\nfactor =', factor)
            print('\npfSORM =', pfsorm)
            print('\nBetaSORM =', betasorm)

        return

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
        matrix = self.Rz[np.ix_(indexes_correlated_xvar, indexes_correlated_xvar)]


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
        L = scipy.linalg.cholesky(self.Rz, lower=True)
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

    def mc(self, nc, ns, delta_lim, nsigma=1.00, igraph=True, iprint=True):
        """
        Monte Carlo Simulation Method
        nc Cycles
        ns Simulations
        Brute force = no adaptive technique

        """
        #
        #
        ti = time.time()
        #
        # Number of variables of the problem
        #
        nc = int(nc)
        ns = int(ns)
        pfc = np.zeros(nc)
        cov_pf = np.zeros(nc)
        pf_mean = np.zeros(nc)
        sum1 = 0.00
        sum2 = 0.00
        fxmax_cycle = np.zeros(nc)
        uk_cycle = np.zeros((ns, self.nxvar))
        #
        # Correlation matrix is self.Rz
        #
        if iprint:
            print('Correlation Matrix after Nataf correction:')
            print(self.Rz)
        #
        # Standard deviation multiplier for MC-IS
        #
        #
        nsigma = 1.00

        #
        #
        # Number of Monte Carlo simulations
        #
        #
        # Matrix xp(ns, self.nxvar) for ns Monte Carlo simulations and self.nxvar random variables
        #
        xp = np.zeros((ns, self.nxvar))
        wp = np.ones(ns)
        fx = np.ones(ns)
        zf = np.zeros((ns, self.nxvar))
        zh = np.zeros((ns, self.nxvar))

        # Matrix dmatrix(ns, self.ndvar) for ns Monte Carlo simulations and self.ndvar design variables

        dmatrix = np.array([self.d.T] * ns)

        #
        # Adaptive cycles
        #

        for icycle in range(nc):
            kcycle = icycle + 1

            #
            # Monte Carlo Simulations
            #
            #
            # Generation of uniform random numbers
            #
            index = icycle % 2
            uk_new = np.random.rand(ns, self.nxvar)
            if index == 0:
                uk_cycle = uk_new.copy()
            else:
                uk_cycle = 1.00 - uk_cycle

            #
            #
            # Step 1 - Generation of the random numbers according to their appropriate distribution
            #

            xp, wp, fx = self.generator(ns)
            
            #
            #
            # Step 2 - Evaluation of the limit state function g(x)
            #
            gx = list(map(self.fel, xp, dmatrix))
            gx = np.array(gx)

            #
            #
            # Step 3 - Evaluation of the indicator function I[g(x)]
            #
            igx = np.where(gx <= 0.00, wp, 0)
            nfail = sum(igx)
            pfc[icycle] = nfail/ns
            sum1 += pfc[icycle]
            sum2 += pfc[icycle] ** 2
            fxmax_cycle[icycle] = fx.max()

            #
            #  Step 6 - Evaluation of the error in the estimation of Pf
            #

            pf_mean[icycle] = sum1 / kcycle
            pf = pf_mean[icycle]
            if pf > 0.00 and kcycle > 1:
                cov_pf[icycle] = 1. / (pf * np.sqrt(kcycle * (kcycle - 1))) * np.sqrt(sum2 - 1. / kcycle * sum1 ** 2)
            else:
                cov_pf[icycle] = 0.00
            delta_pf = cov_pf[icycle]
            # Probability of failure in this cycle
            if iprint: DataVisualize.one_cycle_print_results({'kcycle':kcycle, 'pf':pf, 'delta_pf':delta_pf})
            if delta_pf < delta_lim and kcycle > 3:
                break

        beta = -norm.ppf(pf, 0, 1)
        nsimul = kcycle * ns
        tf = time.time()
        ttotal = tf - ti
        
        # Results viewer   
        if iprint: DataVisualize.print_results({'title':"Monte Carlo - Brute Force", 'beta': beta, 'pf':pf, 'delta_pf': delta_pf, 'nsimul': nsimul, 'gx': gx, 'ttotal': ttotal})

        if igraph: DataVisualize.plot_results({'pf_mean':pf_mean, 'cov_pf':cov_pf, 'kcycle':kcycle})

        return {
            "beta": beta,
            "pf": pf,
            "delta_pf": delta_pf,
            "nsimul": nsimul,
            "ttotal": ttotal
            }
      
    def adaptive(self, nc, ns, delta_lim, nsigma=1.50, igraph=True, iprint=True):
        """
        Monte Carlo Simulations with Importance Sampling (MC-IS)
        Importance sampling with adaptative technique
        Melchers, R.E. Search-based importance sampling.
        Structural Safety, 9 (1990) 117-128

        """
        #
        #
        ti = time.time()
        #
        # Number of variables of the problem
        #
        nfail = 0
        niter = 0
        nc = int(nc)
        ns = int(ns)
        pfc = np.zeros(nc)
        cov_pf = np.zeros(nc)
        pf_mean = np.zeros(nc)
        sum1 = 0.00
        sum2 = 0.00
        fxmax = 0.00
        fxmax_cycle = np.zeros(nc)
        uk_cycle = np.zeros((ns, self.nxvar))

        #
        # Correlation matrix is self.Rz
        #
        if iprint:
            print('Correlation Matrix after Nataf correction:')
            print(self.Rz)

        #
        #
        # Number of Monte Carlo simulations
        #
        #
        # Matrix xp(ns, self.nxvar) for ns Monte Carlo simulations and self.nxvar random variables
        #
        xp = np.zeros((ns, self.nxvar))
        wp = np.ones(ns)
        fx = np.ones(ns)
        zf = np.zeros((ns, self.nxvar))
        zh = np.zeros((ns, self.nxvar))


        # Matrix dmatrix(ns, self.ndvar) for ns Monte Carlo simulations and self.ndvar design variables

        dmatrix = np.array([self.d.T] * ns)

        #
        # Adaptive cycles
        #

        for icycle in range(nc):
            kcycle = icycle + 1

            #
            # Monte Carlo Simulations
            #
            #
            # Step 1 - Generation of the random numbers according to their appropriate distribution
            #

            xp, wp, fx = self.generator(ns)
            #
            #
            # Step 2 - Evaluation of the limit state function g(x)
            #
            gx = list(map(self.fel, xp, dmatrix))
            gx = np.array(gx)

            #
            #
            # Step 3 - Evaluation of the indicator function I[g(x)]
            #
            igx = np.where(gx <= 0.00, wp, 0)
            nfail = sum(igx)
            pfc[icycle] = nfail / ns
            sum1 += pfc[icycle]
            sum2 += pfc[icycle] ** 2
            fxmax_cycle[icycle] = fx.max()

            #
            #  Step 4 - Select adaptative mean
            #
            if nfail == 0:
                #
                # No failures in ns simulations
                #
                imin = np.argmin(gx)
                #
                i = -1
                for var in self.xvar:
                    i += 1
                    var['varhmean'] = xp[imin, i]

            else:
                #
                # Ocurrence of nfail failures in ns simulations
                #
                if fxmax_cycle[icycle] > 1.02 * fxmax:
                    fxmax = fxmax_cycle[icycle]
                    imax = np.argmax(fx)
                    #
                    i = -1
                    for var in self.xvar:
                        i += 1
                        var['varhmean'] = xp[imax, i]

            #
            #  Step 6 - Evaluation of the error in the estimation of Pf
            #

            pf_mean[icycle] = sum1 / kcycle
            pf = pf_mean[icycle]
            if pf > 0.00 and kcycle > 1:
                cov_pf[icycle] = 1. / (pf * np.sqrt(kcycle * (kcycle - 1))) * np.sqrt(sum2 - 1. / kcycle * sum1 ** 2)
            else:
                cov_pf[icycle] = 0.00
            delta_pf = cov_pf[icycle]
            # Probability of failure in this cycle
            if iprint: DataVisualize.one_cycle_print_results({'kcycle':kcycle, 'pf':pf, 'delta_pf':delta_pf})
            if delta_pf < delta_lim and kcycle > 3:
              break

        beta = -norm.ppf(pf, 0, 1)
        nsimul = kcycle * ns
        tf = time.time()
        ttotal = tf - ti
        #
        if iprint: DataVisualize.print_results({'title':"Monte Carlo – Adaptative Importance Sampling", 'beta': beta, 'pf':pf, 'delta_pf': delta_pf, 'nsimul': nsimul, 'gx': gx, 'ttotal': ttotal})

        if igraph: DataVisualize.plot_results({'pf_mean':pf_mean, 'cov_pf':cov_pf, 'kcycle':kcycle})

        return {
            "beta": beta,
            "pf": pf,
            "delta_pf": delta_pf,
            "nsimul": nsimul,
            "ttotal": ttotal
        }

    def bucher(self, nc, ns, delta_lim, nsigma=1.50, igraph=True, iprint=True):
        """
        Monte Carlo Simulations with Importance Sampling (MC-IS)
        Importance sampling with adaptive technique
        BUCHER, C.G. Adaptive sampling – an iterative fast Monte Carlo procedure. Structural
        safety, v. 5, n. 2, p. 119-126, 1988.

        """
        #
        #
        ti = time.time()
        #
        # Number of variables of the problem
        #
        nc = int(nc)
        ns = int(ns)
        xm = np.zeros(self.nxvar)
        sum_xwig = np.zeros(self.nxvar)
        sum_wig = 0.00
        pfc = np.zeros(nc)
        cov_pf = np.zeros(nc)
        pf_mean = np.zeros(nc)
        sum1 = 0.00
        sum2 = 0.00
        uk_cycle = np.zeros((ns, self.nxvar))

        #
        # Correlation matrix is self.Rz
        #
        if iprint:
            print('Correlation Matrix after Nataf correction:')
            print(self.Rz)

        #
        #
        # Number of Monte Carlo simulations
        #
        #
        # Matrix xp(ns, self.nxvar) for ns Monte Carlo simulations and self.nxvar random variables
        #
        xp = np.zeros((ns, self.nxvar))
        wp = np.ones(ns)
        fx = np.ones(ns)

        # Matrix dmatrix(ns, self.ndvar) for ns Monte Carlo simulations and self.ndvar design variables

        dmatrix = np.array([self.d.T] * ns)

        #
        # Adaptive cycles
        #

        for icycle in range(nc):
            kcycle = icycle + 1

            #
            # Monte Carlo Simulations
            #
            # Step 1 - Generation of the random numbers according to their appropriate distribution
            #

            xp, wp, fx = self.generator(ns)
            #
            #
            # Step 2 - Evaluation of the limit state function g(x)
            #
            gx = list(map(self.fel, xp, dmatrix))
            gx = np.array(gx)

            #
            #
            # Step 3 - Evaluation of the indicator function I[g(x)]
            #
            igx = np.where(gx <= 0.00, wp, 0)
            nfail = sum(igx)
            pfc[icycle] = nfail / ns
            sum1 += pfc[icycle]
            sum2 += pfc[icycle] ** 2
            wig = np.copy(igx)

            #
            #  Step 4 - Select adaptive mean
            #
            if nfail == 0:
                #
                # No failures in ns simulations
                #
                imin = np.argmin(gx)
                #
                i = -1
                for var in self.xvar:
                    i += 1
                    xm[i] = xp[imin, i]
                    var['varhmean'] = xm[i]

            else:
                #
                # Ocurrence of nfail failures in ns simulations
                #
                sum_xwig += np.dot(wig.T, xp)
                sum_wig += sum(wig)
                #
                i = -1
                for var in self.xvar:
                    i += 1
                    xm[i] = sum_xwig[i] / sum_wig
                    var['varhmean'] = xm[i]

            #
            #  Step 6 - Evaluation of the error in the estimation of Pf
            #

            pf_mean[icycle] = sum1 / kcycle
            pf = pf_mean[icycle]
            if pf > 0.00 and kcycle > 1:
                cov_pf[icycle] = 1. / (pf * np.sqrt(kcycle * (kcycle - 1))) * np.sqrt(sum2 - 1. / kcycle * sum1 ** 2)
            else:
                cov_pf[icycle] = 0.00
            delta_pf = cov_pf[icycle]
            nc_final = icycle
            # Probability of failure in this cycle
            if iprint: DataVisualize.one_cycle_print_results({'kcycle':kcycle, 'pf':pf, 'delta_pf':delta_pf})
            if delta_pf < delta_lim and kcycle > 3:
              break

        beta = -norm.ppf(pf, 0, 1)
        nsimul = kcycle * ns
        tf = time.time()
        ttotal = tf - ti
        #
        if iprint: DataVisualize.print_results({'title':"Monte Carlo – Bucher Importance Sampling", 'beta': beta, 'pf':pf, 'delta_pf': delta_pf, 'nsimul': nsimul, 'gx': gx, 'ttotal': ttotal})

        if igraph: DataVisualize.plot_results({'pf_mean':pf_mean, 'cov_pf':cov_pf, 'kcycle':kcycle})


        return {
            "beta": beta,
            "pf": pf,
            "delta_pf": delta_pf,
            "nsimul": nsimul,
            "ttotal": ttotal
        }

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
            
            # Correlation matrix is self.Rz
            correlation_matrix = self.Rz

            # Get index of correlated and uncorrelated variables
            index_correlated, index_uncorrelated = self.correlation_summary(correlation_matrix)
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
    
    def sampling_project_point(self, nc, ns, delta_lim, igraph=True, iprint=True):   

      ti = time.time()
      nc = int(nc)
      ns = int(ns)
      pfc = np.zeros(nc)
      sum1 = 0.00
      sum2 = 0.00
      pf_mean = np.zeros(nc)
      cov_pf = np.zeros(nc)

      ## Apriori Results
      formResults = self.form(iHLRF=True, toler=1.e-3, iprint=False)

      ## Project Point
      xk = formResults['xk']
      self.x0 = xk

      ## Using varhmean calculate based on x0
      for var, mean_value in zip(self.xvar, self.x0):
        var['varhmean'] = mean_value

      # Matrix dmatrix(ns, self.ndvar) for ns Monte Carlo simulations and self.ndvar design variables
      dmatrix = np.array([self.d.T] * ns)

      for icycle in range(nc):
        kcycle = icycle + 1

        # Monte Carlo Simulations

        # Generation of uniform random numbers - Antithetic Sampling
        #
        index = icycle % 2
        uk_new = np.random.rand(ns, self.nxvar)
        if index == 0:
            uk_cycle = uk_new.copy()
        else:
            uk_cycle = 1.00 - uk_cycle

        # Step 1 - Generation of the random numbers according to their appropriate distribution
        xp, wp, fx = self.generator(ns)

        # Step 2 - Evaluation of the limit state function g(x)
        gx = list(map(self.fel, xp, dmatrix))
        gx = np.array(gx)

        # Step 3 - Evaluation of the indicator function I[g(x)]
        igx = np.where(gx <= 0.00, wp, 0)
        nfail = sum(igx)
        pfc[icycle] = nfail / ns
        sum1 += pfc[icycle]
        sum2 += pfc[icycle] ** 2
      
        # Step 4 - Evaluation of the error in the estimation of Pf
        pf_mean[icycle] = sum1 / kcycle
        
        pf = pf_mean[icycle]
        if pf > 0.00 and kcycle > 1:
            cov_pf[icycle] = 1. / (pf * np.sqrt(kcycle * (kcycle - 1))) * np.sqrt(sum2 - 1. / kcycle * sum1 ** 2)
        else:
            cov_pf[icycle] = 0.00

        delta_pf = cov_pf[icycle]
        # Plot probability of failure in this cycle
        if iprint: DataVisualize.one_cycle_print_results({'kcycle':kcycle, 'pf':pf, 'delta_pf':delta_pf})
        print('delta_pf',delta_pf)
        if delta_pf < delta_lim and kcycle > 3:
            break

      beta = -norm.ppf(pf, 0, 1)
      nsimul = kcycle * ns
      tf = time.time()
      ttotal = tf - ti

      # Results viewer
      
      if iprint: DataVisualize.print_results({'title':"Monte Carlo – Importance Sampling Based on the Design Point", 'beta': beta, 'pf':pf, 'delta_pf': delta_pf, 'nsimul': nsimul, 'gx': gx, 'ttotal': ttotal})
      if igraph: DataVisualize.plot_results({'pf_mean':pf_mean, 'cov_pf':cov_pf, 'kcycle':kcycle})

      return {
            "beta": beta,
            "pf": pf,
            "delta_pf": delta_pf,
            "nsimul": nsimul,
            "ttotal": ttotal
        }

    def sampling_enhanced(self, nc, ns, delta_lim, igraph=True, iprint=True):

      def calculate_pf(arr_pf, arr_lambda):
            def pf_model(lam, a, b, c, q):
                base = np.maximum(lam - b, 0)
                return q * np.exp(-a * base**c)

            # Initial parameters and limits
            initial_guess = [5.0, 0.1, 1.0, max(arr_pf)]
            bounds = ([0.001, -1.0, 0.1, 0.001], [100.0, 1.5, 5.0, 1.0])

            try:
                params, _ = curve_fit(pf_model, arr_lambda, arr_pf, p0=initial_guess, bounds=bounds, maxfev=20000)
                pf_at_cycle = pf_model(1.0, *params)
                return pf_at_cycle
            
            except Exception as e:
                print("Error in adjustment:", e)
                return np.nan 
                      
      def linspace_between_0_and_1(n):
        """
        Returns a NumPy array with n equally spaced values between 0 and 1,
        excluding 0 and 1.

        Parameters:

        n: int, the number of desired subdivisions in the interval (0, 1)

        Returns:

        numpy.ndarray of floats between 0 and 1, excluding the endpoints.
        """
        return np.linspace(0, 1, n + 2)[1:-1]
      
      ti = time.time()
      nc = int(nc)
      ns = int(ns)

      lambdas = linspace_between_0_and_1(1000)
      nlambdas = len(lambdas)

      pfc = np.zeros((nc, nlambdas))
      pf_mean = np.zeros(nc)
      cov_pf = np.zeros(nc)
      sum1 = 0.00
      sum2 = 0.00

      varhmean_array = [var['varhmean'] for var in self.xvar]
      dvar_array = [var['varvalue'] for var in self.dvar]

      gx_based_varhmean = self.fel(varhmean_array, dvar_array)

      deduction = (1 - lambdas) * gx_based_varhmean
      deductions = np.tile(deduction[:, None], (1, ns))

      # Matrix dmatrix(ns, self.ndvar) for ns Monte Carlo simulations and self.ndvar design variables
      dmatrix = np.array([self.d.T] * ns)

      for icycle in range(nc):
        kcycle = icycle + 1

        # Monte Carlo Simulations
        # Generation of uniform random numbers - Antithetic Sampling
        #
        index = icycle % 2
        uk_new = np.random.rand(ns, self.nxvar)
        if index == 0:
            uk_cycle = uk_new.copy()
        else:
            uk_cycle = 1.00 - uk_cycle

        # Step 1 - Generation of the random numbers according to their appropriate distribution
        xp, wp, fx = self.generator(ns)

        # Step 2 - Evaluation of the limit state function g(x)
        gx = list(map(self.fel, xp, dmatrix))
        gx = np.array(gx)

        # Step 3 - Evaluation of the new limit state function m(x)
        gx_lambdas = np.tile(gx, (nlambdas, 1))
        mx_lambdas = gx_lambdas - deductions

        #Step 4 - Evaluation of the indicator function I[g(x)]
        imx_lambdas = np.where(mx_lambdas <= 0.00, wp, 0)
        nfail = np.sum(imx_lambdas, axis=1)
        pfc[icycle] = nfail / ns       
        pf_cycle = calculate_pf(pfc[icycle], lambdas)
        sum1 += pf_cycle
        sum2 += pf_cycle ** 2

        pf_mean[icycle] = sum1 / kcycle
        pf = pf_mean[icycle]
  
        if pf > 0.00 and kcycle > 1:
            cov_pf[icycle] = 1. / (pf * np.sqrt(kcycle * (kcycle - 1))) * np.sqrt(sum2 - 1. / kcycle * sum1 ** 2)
        else:
            cov_pf[icycle] = 0.00

        delta_pf = cov_pf[icycle]
        
        # Plot probability of failure in this cycle
        if iprint: DataVisualize.one_cycle_print_results({'kcycle':kcycle, 'pf':pf, 'delta_pf':delta_pf})
        if delta_pf < delta_lim and kcycle > 3:
          break

      beta = -norm.ppf(pf, 0, 1)
      nsimul = kcycle * ns
      tf = time.time()
      ttotal = tf - ti

      # Results viewer   
      if iprint: DataVisualize.print_results({'title':"Monte Carlo – Enhanced Importance Sampling", 'beta': beta, 'pf':pf, 'delta_pf': delta_pf, 'nsimul': nsimul, 'gx': gx, 'ttotal': ttotal})

      if igraph: DataVisualize.plot_results({'pf_mean':pf_mean, 'cov_pf':cov_pf, 'kcycle':kcycle})

      return {
            "beta": beta,
            "pf": pf,
            "delta_pf": delta_pf,
            "nsimul": nsimul,
            "ttotal": ttotal
        }
                