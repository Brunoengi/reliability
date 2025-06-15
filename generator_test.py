def var_gen(self, ns, indexes_correlated_xvar, nsigma=1.00, iprint=False):
      """
      Random variables generator for the Monte Carlo Simulation methods
      """

      xvar_correlated = [self.reliability.xvar[i] for i in indexes_correlated_xvar]
      nxvar_correlated = len(xvar_correlated)

      # Gera matrizes e variáveis iniciais
      uk_cycle = np.random.rand(ns, nxvar_correlated)
      x = np.zeros((ns, nxvar_correlated))
      weight = np.ones(ns)
      fxixj = np.ones(ns)
      zf = np.zeros((ns, nxvar_correlated))

      # Submatriz de correlação para variáveis correlacionadas
      matrix = self.reliability.correlation.Rz_rectify[np.ix_(indexes_correlated_xvar, indexes_correlated_xvar)]

      # Cholesky para gerar números gaussianos correlacionados
      L = scipy.linalg.cholesky(matrix, lower=True)
      yk = norm.ppf(uk_cycle)       # Normais independentes
      zk = yk @ L.T                 # Correlacionados

      def update_weights(fx, hx, zf_col, zk_col):
          # Atualiza pesos e produto fxixj para cada variável
          standard_norm_pdf_zf = norm.pdf(zf_col)
          standard_norm_pdf_zk = norm.pdf(zk_col)
          w = (fx / standard_norm_pdf_zf) / (hx / standard_norm_pdf_zk)
          return w, fx / standard_norm_pdf_zf

      for i, var in enumerate(xvar_correlated):
          # Ajusta std se necessário
          if var['varstd'] == 0.0:
              var['varstd'] = float(var['varcov']) * float(var['varmean'])

          namedist = var['vardist'].lower()
          mufx = float(var['varmean'])
          sigmafx = float(var['varstd'])
          muhx = float(var['varhmean'])
          sigmahx = nsigma * sigmafx

          zk_col = zk[:, i]

          if namedist == 'gauss':
              x[:, i] = muhx + sigmahx * zk_col
              fx = norm.pdf(x[:, i], mufx, sigmafx)
              hx = norm.pdf(x[:, i], muhx, sigmahx)
              zf[:, i] = (x[:, i] - mufx) / sigmafx

          elif namedist == 'uniform':
              a = float(var['parameter1'])
              b = float(var['parameter2'])

              def uniform_limits(vars, mux, sigmax):
                  a_, b_ = vars
                  eq1 = (a_ + b_) / 2 - mux
                  eq2 = (b_ - a_) / np.sqrt(12) - sigmax
                  return [eq1, eq2]

              ah, bh = fsolve(uniform_limits, (a, b), args=(muhx, sigmahx))
              uk = norm.cdf(zk_col)
              x[:, i] = ah + (bh - ah) * uk
              zf[:, i] = norm.ppf(uk)
              fx = uniform.pdf(x[:, i], a, b - a)
              hx = uniform.pdf(x[:, i], ah, bh - ah)

          elif namedist == 'lognorm':
              zetafx = sqrt(log(1 + (sigmafx / mufx) ** 2))
              lambdafx = log(mufx) - 0.5 * zetafx ** 2
              zetahx = sqrt(log(1 + (sigmahx / muhx) ** 2))
              lambdahx = log(muhx) - 0.5 * zetahx ** 2
              x[:, i] = np.exp(lambdahx + zk_col * zetahx)
              zf[:, i] = (np.log(x[:, i]) - lambdafx) / zetafx
              fx = norm.pdf(np.log(x[:, i]), lambdafx, zetafx)
              hx = norm.pdf(np.log(x[:, i]), lambdahx, zetahx)

          elif namedist == 'gumbel':
              alphafn = pi / sqrt(6) / sigmafx
              ufn = mufx - euler_gamma / alphafn
              betafn = 1 / alphafn
              alphahn = pi / sqrt(6) / sigmahx
              uhn = muhx - euler_gamma / alphahn
              betahn = 1 / alphahn
              uk = norm.cdf(zk_col)
              x[:, i] = uhn - betahn * np.log(np.log(1 / uk))
              cdfx = gumbel_r.cdf(x[:, i], ufn, betafn)
              zf[:, i] = norm.ppf(cdfx)
              fx = gumbel_r.pdf(x[:, i], ufn, betafn)
              hx = gumbel_r.pdf(x[:, i], uhn, betahn)

          elif namedist == 'frechet':
              deltafx = sigmafx / mufx
              def fkapa(kapa, delta, gsignal):
                  return 1 + delta**2 - gamma(1 + 2 * gsignal / kapa) / (gamma(1 + gsignal / kapa) ** 2)

              kapaf = newton(fkapa, 2.5, args=(deltafx, -1))
              vfn = mufx / gamma(1 - 1 / kapaf)
              deltahx = sigmahx / muhx
              kapah = newton(fkapa, 2.5, args=(deltahx, -1))
              vhn = muhx / gamma(1 - 1 / kapah)
              uk = norm.cdf(zk_col)
              x[:, i] = vhn / (np.log(1 / uk)) ** (1 / kapah)
              ynf = x[:, i] / vfn
              ynh = x[:, i] / vhn
              cdfx = invweibull.cdf(ynf, kapaf)
              zf[:, i] = norm.ppf(cdfx)
              fx = invweibull.pdf(ynf, kapaf) / vfn
              hx = invweibull.pdf(ynh, kapah) / vhn

          elif namedist == 'weibull':
              epsilon = float(var['varinf'])
              deltafx = sigmafx / (mufx - epsilon)
              kapaf = newton(fkapa, 2.5, args=(deltafx, 1))
              w1f = (mufx - epsilon) / gamma(1 + 1 / kapaf) + epsilon
              deltahx = sigmahx / (muhx - epsilon)
              kapah = newton(fkapa, 2.5, args=(deltahx, 1))
              w1h = (muhx - epsilon) / gamma(1 + 1 / kapah) + epsilon
              uk = norm.cdf(zk_col)
              x[:, i] = (w1h - epsilon) * (np.log(1 / (1 - uk))) ** (1 / kapah) + epsilon
              ynf = (x[:, i] - epsilon) / (w1f - epsilon)
              ynh = (x[:, i] - epsilon) / (w1h - epsilon)
              cdfx = weibull_min.cdf(ynf, kapaf)
              zf[:, i] = norm.ppf(cdfx)
              fx = weibull_min.pdf(ynf, kapaf) / (w1f - epsilon)
              hx = weibull_min.pdf(ynh, kapah) / (w1h - epsilon)

          elif namedist == 'beta':
              a = float(var['parameter1'])
              b = float(var['parameter2'])
              q = float(var['parameter3'])
              r = float(var['parameter4'])

              def beta_limits(vars, mux, sigmax, q, r):
                  a_, b_ = vars
                  eq1 = a_ + q / (q + r) * (b_ - a_) - mux
                  eq2 = sqrt((q * r) / (((q + r) ** 2) * (q + r + 1))) * (b_ - a_) - sigmax
                  return [eq1, eq2]

              ah, bh = fsolve(beta_limits, (a, b), args=(muhx, sigmahx, q, r))
              loc, scale = a, b - a
              loch, scaleh = ah, bh - ah
              uk = norm.cdf(zk_col)
              x[:, i] = beta_dist.ppf(uk, q, r, loc=loc, scale=scale)
              fx = beta_dist.pdf(x[:, i], q, r, loc=loc, scale=scale)
              hx = beta_dist.pdf(x[:, i], q, r, loc=loch, scale=scaleh)
              cdfx = beta_dist.cdf(x[:, i], q, r, loc=loc, scale=scale)
              zf[:, i] = norm.ppf(cdfx)

          elif namedist == 'gamma':
              deltafx = sigmafx / mufx
              k = 1 / deltafx ** 2
              v = k / mufx
              a, loc, scale = k, 0, 1 / v

              deltahx = sigmahx / muhx
              kh = 1 / deltahx ** 2
              vh = kh / muhx
              ah, loch, scaleh = kh, 0, 1 / vh

              uk = norm.cdf(zk_col)
              x[:, i] = gamma_dist.ppf(uk, ah, loc=loch, scale=scaleh)
              fx = gamma_dist.pdf(x[:, i], a, loc=loc, scale=scale)
              hx = gamma_dist.pdf(x[:, i], ah, loc=loch, scale=scaleh)
              cdfx = gamma_dist.cdf(x[:, i], a, loc=loc, scale=scale)
              zf[:, i] = norm.ppf(cdfx)

          else:
              raise ValueError(f"Distribuição '{namedist}' não suportada.")

          # Atualiza pesos e produto fxixj para a variável i
          w, fx_over_norm = update_weights(fx, hx, zf[:, i], zk_col)
          weight *= w
          fxixj *= fx_over_norm

      # Correção final com pdf multivariados
      norm_multivarf = multivariate_normal(cov=matrix)
      phif = norm_multivarf.pdf(zf)
      norm_multivarh = multivariate_normal(cov=matrix)
      phih = norm_multivarh.pdf(zk)
      weight *= phif / phih
      fxixj *= phif

      return x, weight, fxixj