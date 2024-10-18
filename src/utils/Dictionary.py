class DictionaryUtils:
  @staticmethod
  def convert_ints_to_floats(dictionary):
    return {key: float(value) if isinstance(value, int) else value for key, value in dictionary.items()}
