from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.base_types.validate_dictionary import ValidateDictionary
from utils.validate.base_types.validate_class import ValidateClass 

class Frechet(AbstractDistribution):
  def __init__(self, props: dict):
    
    ##self.validate_specific_parameters(props)
    super().__init__(props)
    ##ValidateClass.has_invalid_key(self, 'varname', 'vardist', 'varmean', 'varcov', 'varhmean','varstd')

  def validate_specific_parameters(self, props):
    ValidateDictionary.is_dictionary(props)
    ValidateDictionary.has_keys(props, 'varmean', 'varcov')
    ValidateDictionary.check_if_exists(props, 'varcov', lambda d, k: ValidateDictionary.is_greater_or_equal_than(d, k, 0))