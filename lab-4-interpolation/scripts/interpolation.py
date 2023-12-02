#!/usr/bin/python3

# Interpolation of US population using 
# different methods

#======================================
# Imports
#======================================

import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from os import mkdir

#======================================
# Functions
#======================================

#--------------------------------------
# Interpolant in Newtonian form
#--------------------------------------

def build_interpolant_newton(x, y):
  """
  Calculate coeffs of interpolant in Newtonian form.
  """

  coeffs = []
  for k in range(len(x)):
    coeffs.append(div_diff(x, y, k, 0))

  return coeffs

def calculate_newton(x, coeffs, x0) -> float:
  """
  Calculate value of interpolant in Newtonian form in x0.
  """

  value = 0.
  for i in range(len(coeffs)):

    term = coeffs[i]

    for j in range(i):
      term *= (x0 - x[j])

    value += term

  return value

def build_graph_newton(x, coeffs, graph_x):
  """
  Build graph of interpolant in Newtonian form.
  """

  graph_y = [calculate_newton(x, coeffs, point) for point in graph_x]
  build_graph("Interpolant in Newtonian form", "newton", graph_x, graph_y)

def interpolate_newton(x, y, x0: float) -> float:
  
  coeffs = build_interpolant_newton(x, y)
  build_graph_newton(x, coeffs, [i for i in range(int(x[0]), int(x0), 1)])

  return calculate_newton(x, coeffs, x0)

#--------------------------------------
# Spline approximation
#--------------------------------------

def diff(x, y, first, last):
  """
  Helper for calculating right side.
  """ 

  if first == last:
      return y[first]
  return (diff(x, y, first+1, last) - diff(x, y, first, last-1)) / (x[last] - x[first])

def u(x, y, i):
  """
  Helper for calculating right side.
  """  

  return (diff(x, y, i, i + 1) - diff(x, y, i - 1, i)) / (x[i + 1] - x[i - 1])

def f_elem_spline(x, y, order: int, i: int) -> float:
  """
  Returs ith element of right side of the system.
  """

  return 6. * u(x, y, i)

def get_f_spline(x, y, order: int):
  """
  Returs right side of the system.
  """

  return [f_elem_spline(x, y, order, i) for i in range(order)]

def matrix_elem_spline(h: float, x, order: int, i: int, j: int) -> float:
  """
  Returns (i,j) matrix element.
  """

  if i == j:
    return 2.

  if abs(i - j) > 1:
    return 0
  else:
    return 0.5

def get_matrix_spline(h: float, x, y, order: int):
  """
  Returns system's matrix.
  """

  matrix = np.empty((order, order))
  for i in range(order):
    for j in range(order):
      matrix[i][j] = matrix_elem_spline(h, x, order, i, j)

  return matrix

def build_interpolant_spline(h: float, x, y):
  """
  Calculate coeffs of spline approximation.
  """

  order = len(x) - 2

  matrix = get_matrix_spline(h, x, y, order)
  f = get_f_spline(x, y, order)

  return [0, *gauss(matrix, f), 0]

def get_spline_index(x, x0: float) -> int:
  """
  Get index of spline for target point.
  """

  ind = 0

  for i in range(len(x)):
    if x0 <= x[i]:
      ind = i
      break

  return ind

def calculate_spline(h: float, x, y, coeffs, x0) -> float:
  """
  Calculate value of spline approximation in x0.
  """
  
  ind = get_spline_index(x, x0)
  x_i = x[ind]

  a = y[ind]
  b = coeffs[ind] * h / 3 + coeffs[ind-1] * h / 6 + diff(x, y, ind-1, ind)
  c = coeffs[ind]
  d = (coeffs[ind] - coeffs[ind-1]) / h

  return a + b * (x0-x_i) + (c/2) * (x0-x_i)**2 + (d/6) * (x0-x_i)**3

def build_graph_spline(h: float, x, y, coeffs, graph_x):
  """
  Build graph of interpolant in Newtonian form.
  """

  graph_y = [calculate_spline(h, x, y, coeffs, point) for point in graph_x]
  build_graph("Spline approximation", "spline", graph_x, graph_y)

def interpolate_spline(h: float, x, y, x0: float) -> float:

  coeffs = build_interpolant_spline(h, x, y)
  build_graph_spline(h, x, y, coeffs, [i for i in range(int(x[0]), int(x0), 1)])

  return calculate_spline(h, x, y, coeffs, x0)

#--------------------------------------
# Least square method
#--------------------------------------

def f_elem_lsm(x, y, order: int, i: int) -> float:
  """
  Returs ith element of right side of the system.
  """

  power = order - (i + 1)
  return sum([y[k] * pow(x[k], power) for k in range(len(x))])

def get_f_lsm(x, y, order: int):
  """
  Returs right side of the system.
  """

  return [f_elem_lsm(x, y, order, i) for i in range(order)]

def matrix_elem_lsm(x, order: int, i: int, j: int) -> float:
  """
  Returns (i,j) matrix element.
  """

  power = 2 * order - (i + j + 2)
  return sum(pow(x[k], power) for k in range(len(x)))
  
def get_matrix_lsm(x, order: int):
  """
  Returns system's matrix.
  """

  matrix = np.empty((order, order))
  for i in range(order):
    for j in range(order):
      matrix[i][j] = matrix_elem_lsm(x, order, i, j)

  return matrix

def build_interpolant_lsm(x, y, order: int):
  """
  Calculate coeffs of least square method interpolation.
  """

  matrix = get_matrix_lsm(x, order)
  f = get_f_lsm(x, y, order)

  return gauss(matrix, f)

def calculate_lsm(coeffs, x0) -> float:
  """
  Calculate value of least square method interpolant in x0.
  """

  return sum([coeffs[i] * pow(x0, len(coeffs) - (i + 1)) for i in range(len(coeffs))])

def build_graph_lsm(coeffs, graph_x):
  """
  Build graph of least square method interpolant.
  """

  graph_y = [calculate_lsm(coeffs, point) for point in graph_x]
  build_graph("Least square method", "lsm", graph_x, graph_y)

def interpolate_lsm(x, y, x0: float) -> float:

  coeffs = build_interpolant_lsm(x, y, 3)
  build_graph_lsm(coeffs, [i for i in range(int(x[0]), int(x0), 1)])

  return calculate_lsm(coeffs, x0)

#--------------------------------------
# Graphs
#--------------------------------------

def build_graph(title: str, graph_name: str, graph_x, graph_y):

  plt.title(title)

  plt.xlabel("x")
  plt.ylabel("y, interpolation")

  plt.plot(graph_x, graph_y, ".-")

  plt.savefig("res/{}.png".format(graph_name))
  plt.clf()

#--------------------------------------
# SoLAE solver.
#--------------------------------------

def div_diff(x, y, k: int, j: int) -> float:
  """
  Calculate divided difference in x_{j} of order k.
  """

  if k == 0:
    return y[j]

  return (div_diff(x, y, k-1, j+1) - div_diff(x, y, k-1, j)) / (x[j + k] - x[j])

def get_extended_matrix(matrix, f):
  """
  Build extended matrix from matrix and right side of equation
  """

  return np.hstack((matrix,np.array([f]).T))

def gauss(matrix, f):
  """
  Solve SoLAE using Gauss method.
  """

  ext = get_extended_matrix(matrix, f)  

  n = len(matrix)
  for i in range(n):

      leading = i + np.argmax(np.abs(matrix[:,i][i:]))
      ext[[i, leading]] = ext[[leading, i]] 

      ext[i] /= ext[i][i]
      row = ext[i]

      for r in ext[i + 1:]:
          r -= r[i] * row

  for i in range(n - 1, 0, -1):
      row = ext[i]
      for r in reversed(ext[:i]):
          r -= r[i] * row

  return ext[:,-1]

#--------------------------------------
# Misc
#--------------------------------------

def print_res(title: str, x0: float, y_real: float, y0: float):
  
  print(title)
  print(f"In point {x0}: Res == {y0}; Real == {y_real};")
  print(f"Diff == {int(abs(y0 - y_real))};")
  
  print("".join(["#" for i in range(50)]))

#======================================
# Main 
#======================================

def main():

  h = 10.
  x = [1910. + h * i for i in range(10)]
  y = [ 92228496., 106021537., 123202624., 132164569., 151325798., 
       179323175., 203211926., 226545805., 248709873., 281421906. ]
  
  amount = len(y)
  data = {x[i] : y[i] for i in range(amount)}

  x0 = 2010
  y_real = 308745538

  # Interpolant in Newtonian form
  y0 = int(interpolate_newton(x, y, x0))
  print_res("Interpolant in Newtonian form:", x0, y_real, y0)

  # Spline approximation
  y0 = interpolate_spline(h, x, y, x0)
  print_res("Spline approximation:", x0, y_real, y0)

  # Least square method
  y0 = interpolate_lsm(x, y, x0)
  print_res("Least square method:", x0, y_real, y0)

#======================================

if __name__ == '__main__':
  main()
