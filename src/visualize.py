from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

class DataVisualize:
  @staticmethod
  def print_results(inicialTitle, beta, pf, delta_pf, nsimul, gx, ttotal):
      print(f"*** {inicialTitle} ***")
      print(f'\nReliability Index Beta = {beta}')
      print(f'Probability of failure pf = {pf}')
      print(f'COV of pf = {delta_pf}')
      print('nimul = {0:0.4f} '.format(nsimul))
      print(f'Function g(x): mean = {gx.mean()}, std = {gx.std()} ')
      print(f'Processing time = {ttotal} s')

  @staticmethod
  def plot_results(pf_mean, cov_pf, kcycle):
    cycle = np.arange(1, kcycle + 1)

    plt.figure(1, figsize=(8.5, 6))
    plt.plot(cycle, pf_mean[:kcycle])
    plt.title("Convergence of Probability of Failure")
    plt.xlabel("Cycle")
    plt.ylabel("Pf")
    plt.grid(True)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True)) 
    plt.show()

    plt.figure(2, figsize=(8.5, 6))
    plt.plot(cycle, cov_pf[:kcycle])
    plt.title("CoV of the Probability of Failure")
    plt.xlabel("Cycle")
    plt.ylabel("CoV Pf")
    plt.grid(True)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.show()

  @staticmethod
  def one_cycle_print_results(kcycle, pf, delta_pf):
    print('Cycle =', kcycle)
    print(f'Probability of failure pf ={pf}')
    print(f'Coefficient of variation of pf ={delta_pf}')