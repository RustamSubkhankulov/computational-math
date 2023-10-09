#!/usr/bin/python3

#======================================
# Imports
#======================================

import math

import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from os import mkdir

#======================================
# Functions
#======================================

#--------------------------------------
# Генерирование матрицы системы и 
# правой части
#--------------------------------------

def matrix_elem(i: int, j: int):
  """
  Возвращает (i,j) элемент матрицы
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
  Возвращает матрицу системы
  """
  return np.fromfunction(np.vectorize(matrix_elem), (100, 100), dtype=np.double)

def f_elem(i: int):
  """
  Возвращает i-тый элемент правой части системы
  """
  return (i + 1)

def get_f():
  """
  Возвращает правую часть системы
  """
  return np.fromfunction(np.vectorize(f_elem), (100,), dtype = np.double)

def get_extended_matrix(matrix, f):
  """
  Строит расширенную матрицу системы по матрице и правой части
  """
  return np.hstack((matrix,np.array([f]).T))

def norm(vec):
  """
  Норма вектора
  """
  return np.max(np.abs(vec))

def lu_decomposition(matrix):
  """
  Doolittle LU-разложение
  """
  # permutation, lower, upper = sp.linalg.lu(matrix)
  # return lower, upper

  n = len(matrix)

  lower = np.full((n, n), 0.)
  upper = np.full((n, n), 0.)

  for i in range(n):
    lower[i][i] = 1

    if i != 0:
      for j in range(i, n):
        upper[i][j] = matrix[i][j]
        for k in range(0, i):
            upper[i][j] -= lower[i][k] * upper[k][j]

      for j in range(i + 1, n):
        lower[j][i] = matrix[j][i]
        for k in range(0, i):
          lower[j][i] -= lower[j][k] * upper[k][i]
        lower[j][i] /= upper[i][i]
    else:
        for j in range(0, n):
            upper[0, j] = matrix[0, j]
            lower[j, 0] = matrix[j, 0] / matrix[0, 0]

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

def build_graph(title: str, graph_name: str, graph_data):
  """
  Строит и сохраняет график
  """

  plt.title(title)

  plt.xlabel("iteration number")
  
  plt.ylabel("err, log scale")
  plt.yscale("log")

  iterations = [i for i in range(len(graph_data))]
  plt.plot(iterations, graph_data, ".-")

  plt.savefig(f"res/{graph_name}.png")
  plt.clf()

def build_graph_multiple(title: str, graph_name: str, graph_data, labels):
  """
  Строит и сохраняет график с множественным значением серий данных
  """

  plt.title(title)

  plt.xlabel("iteration number")
  
  plt.ylabel("err, log scale")
  plt.yscale("log")

  for ind in range(len(graph_data)):

    iterations = [i for i in range(len(graph_data[ind]))]
    plt.plot(iterations, graph_data[ind], ".-", label=labels[ind])
    plt.legend()

  plt.savefig(f"res/{graph_name}.png")
  plt.clf()

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
def seidel(matrix, f, iterations_number):

  lower, diagonal, upper = ldu_decomposition(matrix)

  B = - np.matmul(np.linalg.inv(lower + diagonal), upper)
  F =   np.matmul(np.linalg.inv(lower + diagonal), f)

  return simple_iteration_method(B, F, matrix, f, iterations_number)

#-----
# Метод Якоби
def jacobi(matrix, f, iterations_number):

  lower, diagonal, upper = ldu_decomposition(matrix)

  B = - np.matmul(np.linalg.inv(diagonal), lower + upper)
  F =   np.matmul(np.linalg.inv(diagonal), f)

  return simple_iteration_method(B, F, matrix, f, iterations_number)

#-----
# Метод Верхней релаксации - Successive Over Relaxation
def successive_over_relaxation(matrix, f, iterations_number, w):

  lower, diagonal, upper = ldu_decomposition(matrix)

  B = - np.matmul(np.linalg.inv(diagonal + w * lower), (w - 1) * diagonal + w * upper)
  F =   np.matmul(np.linalg.inv(diagonal + w * lower), f) * w

  return simple_iteration_method(B, F, matrix, f, iterations_number)

#======================================
# Main 
#======================================

# Точность равенста невязки к нулем
eps = 1e-6

matrix = get_matrix()
f = get_f()

#-----
# Прямые методы
#-----

# Gauss
x = gauss(matrix, f)
assert(norm(np.matmul(matrix,x) - f) < eps)

# LU-decomposition
x = lu(matrix, f)
assert(norm(np.matmul(matrix,x) - f) < eps)

#-----
# Итерационные методы
#-----

res_dir = Path("./res")
if not res_dir.is_dir():
  mkdir(res_dir)

seidel_n = 20
x, err = seidel(matrix, f, seidel_n)
assert(norm(np.matmul(matrix, x) - f) < eps)
build_graph(f"Seidel, number of iterations: {seidel_n}", "seidel", err)

jacobi_n = 45
x, err = jacobi(matrix, f, jacobi_n)
assert(norm(np.matmul(matrix, x) - f) < eps)
build_graph(f"Jacobi, number of iterations: {jacobi_n}", "jacobi", err)

sor_n = 20
errs = []
labels = []

for sor_w in np.arange(1., 2.2, 0.2):

  sor_w = math.ceil(sor_w * 10) / 10
  labels.append("w=" + str(sor_w))

  x, err = successive_over_relaxation(matrix, f, sor_n, sor_w)
  errs.append(err)

build_graph_multiple(
  f"SOR, number of iterations: {sor_n}", 
  f"sor", 
  errs, 
  labels
  )