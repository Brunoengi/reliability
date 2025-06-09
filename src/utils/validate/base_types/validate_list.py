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
  
  @staticmethod
  def uniform_sublength(list):
    try:
      expected_len = len(list[0])
    except TypeError:
      raise TypeError("The first element is not iterable or does not have a length.")

    for i, item in enumerate(list):
      if not hasattr(item, '__len__'):
        raise TypeError(f"Element at index {i} is not iterable or does not have a length.")
      if len(item) != expected_len:
        raise ValueError(
          f"Element at index {i} has length {len(item)}, "
          f"expected length {expected_len}."
          )