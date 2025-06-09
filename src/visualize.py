from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

class DataVisualize:
  @staticmethod
  def print_results(result: dict):
    # Print the title if available
    title = result.get('title')
    if title:
      print(f"*** {title} ***")

    # Dispatch table: maps keys to corresponding print functions
    printers = {
      'beta': lambda v: print(f'\nReliability Index (β) = {v}'),
      'pf': lambda v: print(f'Probability of Failure (Pf) = {v}'),
      'delta_pf': lambda v: print(f'Coefficient of Variation of Pf = {v}'),
      'nsimul': lambda v: print(f'Number of Simulations = {v:0.4f}'),
      'gx': lambda v: print(f'Limit State Function g(x): mean = {v.mean()}, std = {v.std()}') 
          if hasattr(v, 'mean') and hasattr(v, 'std') 
          else print('Limit State Function g(x): invalid data'),
      'ttotal': lambda v: print(f'Processing Time = {v} s')
    }

    # Execute each print function if the corresponding key is present
    for key, printer in printers.items():
      if key in result:
        try:
            printer(result[key])
        except Exception as e:
            print(f"[Error printing '{key}']: {e}")

  @staticmethod
  def plot_results(results: dict):
    pf_mean = results.get("pf_mean")
    cov_pf = results.get("cov_pf")
    kcycle = results.get("kcycle")

    if pf_mean is None or cov_pf is None or kcycle is None:
      print("Missing required data: 'pf_mean', 'cov_pf', or 'kcycle'.")
      return

    cycle = np.arange(1, kcycle + 1)

    # Plot Probability of Failure convergence
    plt.figure(1, figsize=(8.5, 6))
    plt.plot(cycle, pf_mean[:kcycle])
    plt.title("Convergence of the Probability of Failure")
    plt.xlabel("Cycle")
    plt.ylabel("Pf")
    plt.grid(True)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.show()

    # Plot Coefficient of Variation of Pf
    plt.figure(2, figsize=(8.5, 6))
    plt.plot(cycle, cov_pf[:kcycle])
    plt.title("Coefficient of Variation of Pf")
    plt.xlabel("Cycle")
    plt.ylabel("CoV of Pf")
    plt.grid(True)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.show()

  @staticmethod
  def one_cycle_print_results(result: dict):
    # Dispatch table
    printers = {
        'kcycle': lambda v: print(f'Cycle = {v}'),
        'pf': lambda v: print(f'Probability of Failure (Pf) = {v}'),
        'delta_pf': lambda v: print(f'Coefficient of Variation of Pf = {v}')
      }

    # Execute print functions for available keys
    for key, printer in printers.items():
      if key in result:
        try:
          printer(result[key])
        except Exception as e:
          print(f"[Error printing '{key}']: {e}")