from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.base_types.validate_dictionary import ValidateDictionary
from utils.validate.base_types.validate_class import ValidateClass 

class LogNormal(AbstractDistribution):
  def __init__(self, dictionaryInfo):

    ValidateDictionary.is_dictionary(dictionaryInfo)
    ValidateDictionary.check_possible_arrays_keys(dictionaryInfo,['varmean','varstd'], ['varmean','varcov'])
    super().__init__(dictionaryInfo)
    ValidateClass.has_invalid_key(self, 'varname', 'vardist', 'varmean', 'varstd', 'varcov', 'varhmean')