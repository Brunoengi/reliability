import inspect

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