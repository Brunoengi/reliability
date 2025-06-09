from utils.validate.base_types.validate_dictionary import ValidateDictionary 
from utils.validate.base_types.validate_list import ValidateList

class ValidateDvar():
  def __init__(self, dvar):
    ValidateList.is_list(dvar)
    ValidateList.same_type_elements(dvar)

    for item in dvar:
      ValidateDictionary.is_dictionary(item)
      ValidateDictionary.has_keys(item, 'varname', 'varvalue')
      ValidateDictionary.has_invalid_keys(item, 'varname', 'varvalue')