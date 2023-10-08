#!/usr/bin/python3

#======================================
# Imports
#======================================

import math

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

#======================================
# Functions
#======================================

#--------------------------------------
# Генерирование матрицы системы и 
# правой части
#--------------------------------------

def matrix_elem(i: int, j: int):
  """
  Returns a(i,j) - coeff in i'th row and j'th column
  """
  if i == 99:
    return 1
  else:
    if i == j:
      return 10
    elif i == (j + 1) or j == (i + 1):
      return 1
    else: 
      return 0

def get_matrix():
  """
  Returns matrix of the the system
  """
  return np.fromfunction(np.vectorize(matrix_elem), (100, 100), dtype=np.double)

def f_elem(i: int):
  """
  Returns f(i) - coeff on right side of system in i'th row
  """
  return (i + 1)

def get_f():
  """
  Returns right side of the system
  """
  return np.fromfunction(np.vectorize(f_elem), (100,), dtype = np.double)

def get_extended_matrix(matrix, f):
  """
  Returns extended matrix of the system
  """
  return np.hstack((matrix,np.array([f]).T))

def norm(vec):
  """
  Норма вектора
  """
  return np.max(np.abs(vec))

def lu_decomposition(matrix):
 
  n = len(matrix)

  lower = [[0 for x in range(n)] for y in range(n)]
  upper = [[0 for x in range(n)] for y in range(n)]

  for i in range(n):
    for k in range(i, n):
      sum = 0
      for j in range(i):
          sum += (lower[i][j] * upper[j][k])
      upper[i][k] = matrix[i][k] - sum

    for k in range(i, n):
      if (i == k):
          lower[i][i] = 1
      else:
        sum = 0
        for j in range(i):
            sum += (lower[k][j] * upper[j][i])

        lower[k][i] = int((matrix[k][i] - sum) / upper[i][i])

  return lower, upper

#--------------------------------------
# Прямые методы
#--------------------------------------

def gauss(matrix, f):

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

def lu(matrix, f):

  n = len(matrix)
  lower, upper = lu_decomposition(matrix)

  y = [0 for i in range(n)]
  for i in range(n):
    y[i] = f[i] - np.matmul(lower[i], y)

  x = [0 for i in range(n)]
  for i in range(n - 1, -1, -1):
    x[i] = (y[i]) - np.matmul(upper[i], x) / upper[i][i]

  return x;

#======================================
# Main 
#======================================

eps = 10**(-6)

matrix = get_matrix()
f = get_f()

#-----
# Gauss
x = gauss(matrix, f)
assert norm(np.matmul(matrix,x) - f) < eps 

#-----
# LU-decomposition
x = lu(matrix, f)
assert norm(np.matmul(matrix,x) - f) < eps