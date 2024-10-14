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
    def is_float(dictionary, key):
      if not isinstance(dictionary[key], float):
          raise TypeError(f"Error: The key '{key}' must match with a float type value")
            
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
    def check_keys_count(dictionary, n, *possibleKeys):
        count = sum(1 for key in possibleKeys if key in dictionary)
        if count != n:
            raise ValueError(f"Error: Exactly {n} of the following keys must be present: {possibleKeys}. Found {count}.")
     
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
                f"Sets present in the dictionary: {valid_str or 'Nenhum'}"
            )
