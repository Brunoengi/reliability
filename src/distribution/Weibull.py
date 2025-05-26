from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.Dictionary import ValidateDictionary
from utils.validate.Class import ValidateClass

class Weibull(AbstractDistribution):
  def __init__(self, dictionaryInfo):

    ValidateDictionary.is_dictionary(dictionaryInfo)
    ValidateDictionary.has_keys(dictionaryInfo, 'varmean', 'varstd', 'parameter3')
    self.varinf = float(dictionaryInfo['parameter3'])
    ValidateClass.has_invalid_key(self, 'varname', 'vardist', 'varmean', 'varstd', 'parameter3', 'varhmean','varinf')
