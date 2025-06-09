from utils.validate.base_types.validate_function import ValidateFunction


class ValidateGx():

  def __init__(self, gx):

    ValidateFunction.is_function(gx)
    ValidateFunction.has_n_args(gx, 2)
    ValidateFunction.args_are_lists(gx, 2)

 
  