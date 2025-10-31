from ..base_types.validate_function import ValidateFunction
import numpy as np


class ValidateGx():

  def __init__(self, gx, nxvar: int, ndvar: int):

    ValidateFunction.is_function(gx)
    ValidateFunction.has_at_least_n_args(gx, 2)
    ValidateFunction.first_n_args_are_lists(gx, 2)
    self.validate_return(gx, nxvar, ndvar)

  def validate_return(self, gx, nxvar, ndvar):
    try:
      x = [1.0] * nxvar
      d = [1.0] * ndvar
      result = gx(x, d)

      if not isinstance(result, (int, float)):
        raise TypeError(
          f"The function should return a numeric value, got {type(result).__name__}."
        )

    except IndexError as e:
      raise RuntimeError(
        f"The function raised an IndexError, likely due to accessing x[i] or d[i] "
        f"with lengths nxvar={nxvar} and ndvar={ndvar}: {e}"
      )
    except Exception as e:
      raise RuntimeError(f"The function raised an error when called: {e}")