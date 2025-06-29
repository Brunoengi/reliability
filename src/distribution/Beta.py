from distribution.AbstractDistribution import AbstractDistribution 
from utils.validate.base_types.validate_dictionary import ValidateDictionary
from utils.validate.base_types.validate_class import ValidateClass 
import numpy as np

class Beta(AbstractDistribution):
  def __init__(self, props):

    self.validate_specific_parameters(props)

    a = props['parameter1']
    b = props['parameter2']
    q = props['parameter3']
    r = props['parameter4']

    self.varmean = float(a + q / (q + r) * (b - a))
    self.varstd = float(np.sqrt((q * r) / ((q + r) **2 * (q + r + 1)) * (b - a) ** 2))
    self.varhmean = float(self.varmean)

    super().__init__(props)
    ValidateClass.has_invalid_key(self, 'varname', 'vardist','varmean','varcov','varstd','varhmean','parameter1','parameter2','parameter3','parameter4')

  def validate_specific_parameters(self, props):
    ValidateDictionary.is_dictionary(props)
    ValidateDictionary.check_possible_arrays_keys(props, ['parameter1', 'parameter2', 'parameter3', 'parameter4'])