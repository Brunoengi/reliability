from ..base_types.validate_dictionary import ValidateDictionary

# Contains general validations, specific validations are inserted in each distribution class

class ValidateXvar():
  def __init__(self, xvar: dict):
    ValidateDictionary.has_keys(xvar, 'varname','vardist')
    ValidateDictionary.is_string(xvar, 'varname', 'vardist')
