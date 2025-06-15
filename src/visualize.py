from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

class DataVisualizeSimulation:
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

class DataVisualizeTransformation:
  @staticmethod
  def one_cycle_print_results(result: dict):
      print("\n===== TRANSFORMATION RESULTS =====")

      # Scalar variables
      scalar_printers = {
          'kiter': lambda v: print(f'Iteration number = {v:d}'),
          'gxk': lambda v: print(f'g(x) = {v:0.5e}'),
          'erro1': lambda v: print(f'Error in direction (erro1) = {v:0.5e}'),
          'erro2': lambda v: print(f'Error in limit state (erro2) = {v:0.5e}'),
          'beta': lambda v: print(f'Reliability Index (Beta) = {v:0.4f}'),
          'ck': lambda v: print(f'Step size coefficient (ck) = {v:0.5e}'),
          'lambdak': lambda v: print(f'Line search factor (lambda_k) = {v:0.5e}')
      }

      for key, printer in scalar_printers.items():
          if key in result and np.isfinite(result[key]):
              try:
                  printer(result[key])
              except Exception as e:
                  print(f"[Error printing '{key}']: {e}")

      # DataFrame of variable values
      df_keys = ['namevar', 'dist', 'muxneqk', 'sigmaxneqk', 'xk', 'yk', 'alpha']
      if all(k in result for k in df_keys):
          try:
              datadict = {
                  'xvar': result['namevar'],
                  'prob_dist': result['dist'],
                  'mean': result['muxneqk'],
                  'std_dev': result['sigmaxneqk'],
                  'xk': result['xk'],
                  'yk': result['yk'],
                  'alpha': result['alpha']
              }
              data = pd.DataFrame(datadict)
              print("\n--- Variables Table ---")
              print(data)
          except Exception as e:
              print(f"[Error creating DataFrame]: {e}")
      else:
          print("\n[Info]: Insufficient data to create the design variable table.")

      print("===================================")

  @staticmethod
  def print_results(result: dict):
        """
        Prints final reliability results from FORM/SORM transformations,
        such as reliability index (Beta) and probability of failure (Pf).
        """
        printers = {
            'pf': lambda v: print(f'\nProbability of Failure Pf = {v:0.4e}'),
            'beta': lambda v: print(f'\nReliability Index β = {v:0.4e}')
        }

        for key, printer in printers.items():
            if key in result:
                try:
                    printer(result[key])
                except Exception as e:
                    print(f"[Error printing '{key}']: {e}")