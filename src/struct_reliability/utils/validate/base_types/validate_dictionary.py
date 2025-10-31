# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 19:41:47 2024

@author: BrunoTeixeira
"""

class ValidateDictionary:
  @staticmethod
  def has_keys(dictionary, *keys):
    try:
      return all(key in dictionary for key in keys)
    except KeyError:
      missing_keys = [key for key in keys if key not in dictionary]
      raise KeyError(f"Error: The key(s) '{missing_keys}' are missing in {dictionary}")  
  
  @staticmethod
  def is_float(dictionary, *keys):
    for key in keys:
      value = dictionary.get(key)
      if not isinstance(value, float):
        raise TypeError(f"Error: The key '{key}' must be associated with a float type value, but got {type(value).__name__} with value {value}")
          
  @staticmethod
  def is_dictionary(dictionary):
    try:
      return isinstance(dictionary, dict)
    except TypeError:
      raise TypeError("Error: You need to insert a dictionary")
  
  @staticmethod
  def is_greater_than(dictionary, key, value): 
    if(float(dictionary[key]) <= value):
        raise ValueError(f"Error: The value of '{key}' must be greater than {value}. Current value is {dictionary[key]}.")
    
  @staticmethod
  def is_greater_or_equal_than(dictionary, key, value): 
    if(float(dictionary[key]) < value):
        raise ValueError(f"Error: The value of '{key}' must be greater or equal than {value}. Current value is {dictionary[key]}.")
      
  @staticmethod
  def check_keys_count(dictionary, n, *possible_keys):
    count = sum(1 for key in possible_keys if key in dictionary)
    if count != n:
        raise ValueError(f"Error: Exactly {n} of the following keys must be present: {possible_keys}. Found {count}.")
    
  @staticmethod
  def check_possible_arrays_keys(dictionary, *arrays):
    dictionary_keys = set(dictionary.keys())
    valid_arrays = [array for array in arrays if dictionary_keys.issuperset(set(array))]
    if len(valid_arrays) != 1:
      array_str = ', '.join(str(array) for array in arrays)
      valid_str = ', '.join(str(array) for array in valid_arrays)
      raise ValueError(
        f"Error: Exactly one of the key sets must be present in the dictionary.\n"
        f"Sets supplied: {array_str}\n"
        f"Sets present in the dictionary: {valid_str or 'None'}"
      )
    
  @staticmethod
  def has_invalid_keys(dictionary, *all_valid_keys):
    if not isinstance(dictionary, dict):
      raise TypeError(f"Expected a dictionary but received '{type(dictionary).__name__}'.")

    invalid_keys = [key for key in dictionary if key not in all_valid_keys]

    if invalid_keys:
      raise ValueError(
        f"Invalid keys found: {invalid_keys}. "
        f"Valid keys are: {list(all_valid_keys)}."
      )

  @staticmethod
  def is_string(dictionary, *keys):
    for key in keys:
      value = dictionary.get(key)
      if not isinstance(value, str):
        raise TypeError(f"Error: The key '{key}' must be associated with a string type value, but got {type(value).__name__} with value {value}")

  @staticmethod 
  def check_if_exists(dictionary, key, check_function):
    if key in dictionary:
      check_function(dictionary, key)