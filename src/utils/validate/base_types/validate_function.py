import inspect
import numpy as np

class ValidateFunction:
  
  @staticmethod
  def is_function(function):
    if not callable(function):
      raise TypeError("The given argument is not a valid function.")
    
  @staticmethod
  def has_n_args(function, expected_arg_count):
    # Get the function's parameters
    sig = inspect.signature(function)
    params = sig.parameters.values()

    # Filter only fixed arguments (ignore *args and **kwargs)
    fixed_params = [
      p for p in params 
      if p.kind in (
        inspect.Parameter.POSITIONAL_ONLY,
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
        inspect.Parameter.KEYWORD_ONLY
      )
    ]

    if len(fixed_params) != expected_arg_count:
      raise ValueError(
      f"The function must have exactly {expected_arg_count} explicit argument(s). "
      f"Found: {[p.name for p in fixed_params]}"
      )

  @staticmethod
  def has_at_least_n_args(function, min_arg_count):
    # Get the function's parameters
    sig = inspect.signature(function)
    params = sig.parameters.values()

    # Filter only fixed arguments (ignore *args and **kwargs)
    fixed_params = [
      p for p in params 
      if p.kind in (
        inspect.Parameter.POSITIONAL_ONLY,
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
        inspect.Parameter.KEYWORD_ONLY
        )
    ]

    if len(fixed_params) < min_arg_count:
      raise ValueError(
        f"The function must have at least {min_arg_count} explicit argument(s). "
        f"Found: {[p.name for p in fixed_params]}"
      )

  @staticmethod
  def args_are_lists(function, expected_arg_count):
    # Create test lists (simple lists)
    test_args = [[0, 1, 2, 3] for _ in range(expected_arg_count)]

    try:
      # Try to execute the function with the test lists
      function(*test_args)
    except Exception as e:
      raise TypeError(
          f"The function raised an error when called with list arguments: {e}"
      )

  @staticmethod
  def first_n_args_are_lists(function, n_required_lists):
    # Check if the input is a callable function
    if not callable(function):
        raise TypeError("The given object is not a callable function.")

    # Get the function's signature and parameters
    sig = inspect.signature(function)
    params = list(sig.parameters.values())
    total_params = len(params)

    # Ensure the function has at least n_required_lists parameters
    if n_required_lists > total_params:
      raise ValueError(
        f"The function only accepts {total_params} argument(s), "
        f"but {n_required_lists} list argument(s) were required."
      )

    # Prepare test arguments:
    # Pass empty lists to the first n_required_lists parameters to check type acceptance
    test_args = [[] for _ in range(n_required_lists)]

    # Fill remaining parameters with None (or other dummy value)
    for _ in range(total_params - n_required_lists):
      test_args.append(None)

    try:
      # Call the function with the test arguments
      function(*test_args)
    except IndexError:
      # Ignore IndexError, assuming it's caused by empty lists being too short
      # This means the function accepted lists as input but failed due to length
      pass
    except Exception as e:
      # For any other exceptions, raise a TypeError with details
      raise TypeError(
        f"The function raised an error when called with list arguments "
        f"for the first {n_required_lists} parameter(s): {e}"
      )

  @staticmethod
  def gfunction(x, d):
    # Validate input types
    if not isinstance(x, (list, tuple, np.ndarray)):
      raise TypeError(f"'x' must be a list, tuple, or numpy.ndarray, but got {type(x)}.")
    if not isinstance(d, (list, tuple, np.ndarray)):
      raise TypeError(f"'d' must be a list, tuple, or numpy.ndarray, but got {type(d)}.")
