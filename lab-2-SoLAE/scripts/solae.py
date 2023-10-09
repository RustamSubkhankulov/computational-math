#!/usr/bin/python3

#======================================
# Imports
#======================================

import math
import scipy

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
  """
  LU-разложение
  """
  permutation, lower, upper = scipy.linalg.lu(matrix)
  return lower, upper

def ldu_decomposition(matrix):
  """
  LDU-разложение  
  """
  u = np.triu(matrix)
  l = np.tril(matrix)

  lower = matrix - u
  diagonal = l + u - matrix
  upper = matrix - l

  return lower, diagonal, upper

#--------------------------------------
# Прямые методы
#--------------------------------------

#-----
# Метод Гаусса
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

#-----
# LU-разложение
def lu(matrix, f):

  n = len(matrix)
  lower, upper = lu_decomposition(matrix)

  x = [0 for i in range(n)]
  y = [0 for i in range(n)]

  for i in range(1, n + 1):

    summ = 0
    for k in range(1, i):
      summ += lower[i-1][k-1] * y[k-1]

    y[i-1] = f[i-1] - summ

  for i in range(n, 0, -1):

    summ = 0
    for k in range(n, i, -1):
      summ += upper[i-1][k-1] * x[k-1]

    x[i-1] = (y[i-1] - summ) / upper[i-1][i-1]

  return x;

#--------------------------------------
# Итерационные методы
#--------------------------------------

#-----
# Метод простой итерации
def simple_iteration_method(B, F, matrix, f, iterations_number):

  n = len(F)

  x = [0 for i in range(n)]
  err = [0 for i in range(iterations_number)]

  for iter in range(iterations_number):
    
    x = np.matmul(B, x) + F
    
    diff =  f - np.matmul(matrix, x)
    err[iter] = norm(diff) 

  return x, err

#-----
# Метод Зейделя
def seidel(matrix, f):

  lower, diagonal, upper = ldu_decomposition(matrix)

  B = -np.matmul(np.linalg.inv(lower + diagonal), upper)
  F =  np.matmul(np.linalg.inv(lower + diagonal), f)

#-----
# Метод Якоби
# def jacobi(matrix, f):


#-----
# Метод Зейделя
# def upper_relaxation(matrix, f):

#======================================
# Main 
#======================================

eps = 1e-6

matrix = get_matrix()
f = get_f()

#-----
# Прямые методы
#-----

#-----
# Gauss
x = gauss(matrix, f)
assert norm(np.matmul(matrix,x) - f) < eps 

#-----
# LU-decomposition
x = lu(matrix, f)
assert norm(np.matmul(matrix,x) - f) < eps

#-----
# Итерационные методы методы
#-----
