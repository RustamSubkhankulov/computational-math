#!/usr/bin/python3

# System of non-linear equations
# sin(x+2) - y   = 1.5
# x + cos(y - 2) = 0.5

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
# System of non-linear equations
#--------------------------------------

# sin(x+2) - y   = 1.5
# x + cos(y - 2) = 0.5

#--------------------------------------
# Simple iteration method
#--------------------------------------

def x_next_iteration_sim(x_n: float, y_n: float) -> float:
  """
  x(n+1) = 0.5 - cos(y(n) - 2)
  """
  return 0.5 - math.cos(y_n - 2.0)

def y_next_iteration_sim(x_n: float, y_n: float) -> float:
  """
  y(n+1) = sin(x(n) + 2) - 1.5
  """
  return math.sin(x_n + 2.0) - 1.5

#--------------------------------------
# Newton's method
#--------------------------------------

def F(x: float, y: float):
  """
  [sin(x+2) - y - 1.5; x + cos(y - 2) - 0.5]
  """
  return [math.sin(x + 2.0) - y - 1.5, x + math.cos(y - 2.0) - 0.5]

def J(x: float, y: float):
  """
  Jacoby matrix of the system
  """
  j = np.zeros((2,2))

  j[0][0] = +math.cos(x + 2.0)
  j[0][1] = -1.0
  j[1][0] = +1.0
  j[1][1] = -math.sin(y - 2.0)

  return j

def x_next_iteration_newton(x_n: float, y_n: float) -> float:

  j = J(x_n, y_n)
  return x_n - (np.linalg.inv(j) @ F(x_n, y_n))[0]

def y_next_iteration_newton(x_n: float, y_n: float) -> float:
  
  j = J(x_n, y_n)
  return y_n - (np.linalg.inv(j) @ F(x_n, y_n))[1]

#--------------------------------------
# Iterate
#--------------------------------------

def iterate(x_0: float, y_0: float, x_next_iteration, y_next_iteration, eps: float):
  """
  Make iterations until desired presicion is achieved.
  """

  iterations_x = [x_0]
  iterations_y = [y_0]

  prev_x = x_0
  prev_y = y_0

  next_x = x_next_iteration(prev_x, prev_y)
  next_y = y_next_iteration(prev_x, prev_y)

  iterations_x.append(next_x)
  iterations_y.append(next_y)

  while abs(next_x - prev_x) > eps or abs(next_y - prev_y) > eps:

    prev_x = next_x
    prev_y = next_y

    next_x = x_next_iteration(prev_x, prev_y)
    next_y = y_next_iteration(prev_x, prev_y)

    iterations_x.append(next_x)
    iterations_y.append(next_y)

  return (iterations_x, iterations_y)

#--------------------------------------
# Graph build
#--------------------------------------

def build_graph(title: str, graph_name: str, iterations_x, iterations_y):
  """
  Build and generate graph as img.
  """
  plt.title(title)

  plt.xlabel("n, iteration number")
  plt.ylabel("n-th iteration")

  numbers_x = [i for i in range(len(iterations_x))]
  plt.plot(numbers_x, iterations_x, ".-", label="x(n)")

  numbers_y = [i for i in range(len(iterations_y))]
  plt.plot(numbers_y, iterations_y, ".-", label="y(n)")

  plt.savefig(f"res/{graph_name}.png")
  plt.clf()

#======================================
# Main 
#======================================

# Precision
eps = 1e-6

# Real root values
x_real = +1.346
y_real = -1,703

#--------------------------------------
# Simple iteration method
#--------------------------------------

# start approximations
x_start = +1.0  
y_start = -2.0  

(iterations_x, iterations_y) = iterate(
  x_start, 
  y_start, 
  x_next_iteration_sim, 
  y_next_iteration_sim, 
  eps
  )

build_graph(
  f"Simple iteration method: x=~{x_real}, y=~{y_real}", 
  "SIM_SYS", 
  iterations_x, 
  iterations_y
  )

#--------------------------------------
# Newton's method
#--------------------------------------

# start approximations
x_start = +1.0  
y_start = -2.0  

(iterations_x, iterations_y) = iterate(
  x_start, 
  y_start, 
  x_next_iteration_newton, 
  y_next_iteration_newton, 
  eps
  )

build_graph(
  f"Newton's method: x=~{x_real}, y=~{y_real}", 
  "NEWTON_SYS", 
  iterations_x, 
  iterations_y
  )