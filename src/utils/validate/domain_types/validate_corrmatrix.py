from utils.validate.base_types.validate_list import ValidateList
import math
import numpy as np


class ValidateCorrelationMatrix():
  def __init__(self, matrix):
    self.matrix = matrix

    self.is_valid()
    ValidateList.same_type_elements(matrix)
    self.validate_diagonal_is_ones()
    self.validate_symmetry()

  def validate_diagonal_is_ones(self):
    size = len(self.matrix)
    for i in range(size):
      try:
        value = self.matrix[i][i]
      except IndexError:
        raise ValueError(f"Matrix is not square: row {i} does not have enough elements.")

      if value != 1.0:
        raise ValueError(
          f"Invalid correlation matrix: diagonal element at position ({i},{i}) is {value}, expected 1.0."
        )
      
  def validate_symmetry(self):
    size = len(self.matrix)
    for i in range(size):
      for j in range(i + 1, size):
        a = self.matrix[i][j]
        b = self.matrix[j][i]
        if not math.isclose(a, b, rel_tol=1e-9):
          raise ValueError(
            f"Matrix is not symmetric: element at ({i},{j}) = {a} "
            f"differs from ({j},{i}) = {b}."
          )
  
  def is_valid(self):     
    if not isinstance(self.matrix, (list, np.ndarray)):
      raise TypeError(f"Expected an object of type 'list' or 'ndarray', but received '{type(self.matrix).__name__}'.")