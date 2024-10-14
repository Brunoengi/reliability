from src.distribution.AbstractDistribution import AbstractDistribution 
from src.utils.validate.Dictionary import ValidateDictionary
from src.utils.validate.Class import ValidateClass

class Weibull(AbstractDistribution):
  def __init__(self, dictionaryInfo):

    ValidateDictionary.is_dictionary(dictionaryInfo)
    ValidateDictionary.has_keys(dictionaryInfo, 'varmean', 'varstd', 'parameter3')
    self.varinf = float(dictionaryInfo['parameter3'])
    ValidateClass(self, 'varname', 'vardist', 'varmean', 'varstd', 'parameter3', 'varhmean')
    