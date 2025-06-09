class ValidateList():
  @staticmethod
  def is_list(obj):
    if not isinstance(obj, list):
      raise TypeError(f"Expected an object of type 'list', but received '{type(obj).__name__}'.")
  
  @staticmethod
  def same_type_elements(list):
    first_type = type(list[0])
    for i, elem in enumerate(list):
      if not isinstance(elem, first_type):
        raise ValueError(
          f"Element in index {i} is the type '{type(elem).__name__}', "
          f"expected '{first_type.__name__}'."
        )