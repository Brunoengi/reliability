from src.utils.Exceptions import InvalidKeyError

class ValidateClass:
  @staticmethod
  def has_invalid_key(the_class, *all_valid_keys):
      # Gets the object's keys (the instance attributes)
      actual_keys = set(the_class.__dict__.keys())
      
      # Converts all_valid_keys to a set for comparison
      valid_keys = set(all_valid_keys)

      # Checks for invalid keys
      invalid_keys = actual_keys - valid_keys
      
      if invalid_keys:
          raise InvalidKeyError(invalid_keys)