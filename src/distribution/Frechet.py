from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.Dictionary import ValidateDictionary
from utils.validate.Class import ValidateClass

class Frechet(AbstractDistribution):
  def __init(self, dictionaryInfo):

    ValidateDictionary.is_dictionary(dictionaryInfo)
    ValidateDictionary.has_keys(dictionaryInfo, 'varmean', 'varcov')
    super().__init__(dictionaryInfo)
    ValidateClass.has_invalid_key(self, 'varname', 'vardist', 'varmean', 'varcov', 'varhmean')