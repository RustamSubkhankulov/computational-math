#!/usr/bin/python3

#======================================
# Imports
#======================================

import math

import matplotlib.pyplot as plt

from pathlib import Path
from os import mkdir

#======================================
# Functions
#======================================

#--------------------------------------
# Non-linear equation
#--------------------------------------

def f(x: float) -> float:
  """
  2 * x ^2 + 5 * x - 3
  """
  return 2.0 * pow(x, 2.0) + 5.0 * x - 3.0

def f_deriv(x: float) -> float:
  """
  4 * x + 5
  """
  return 4.0 * x + 5.0

#--------------------------------------
# Simple iteration method
#--------------------------------------

def x_next_iteration_sim(x_n: float) -> float:
  """
  x(n+1) = 0,6 - 0,4 * (x(n))^2
  """
  return 0.6 - 0.4 * pow(x_n, 2.0)

#--------------------------------------
# Newton's method
#--------------------------------------

def x_next_iteration_newton(x_n: float) -> float:
  """
  x(n+1) = x(n) - f(x(n)) / f'(x(n))
  """
  return x_n - (f(x_n)) / (f_deriv(x_n))

#--------------------------------------
# Iterate
#--------------------------------------

def iterate(x_0: float, x_next_iteration, eps: float):
  """
  Make iterations until desired presicion is achieved.
  """

  iterations = [x_0]

  prev = x_0

  next = x_next_iteration(prev)
  print(next)
  iterations.append(next)

  while abs(next - prev) > eps:

    prev = next 

    next = x_next_iteration(prev)
    print(next)
    iterations.append(next)

  return iterations

#--------------------------------------
# Graph build
#--------------------------------------

def build_graph(title: str, graph_name: str, iterations):
  """
  Build and generate graph as img.
  """
  plt.title(title)

  plt.xlabel("n, iteration number")
  plt.ylabel("x(n), n-th iteration")

  numbers = [i for i in range(len(iterations))]
  plt.plot(numbers, iterations, ".-")

  plt.savefig(f"res/{graph_name}.png")
  plt.clf()

#======================================
# Main 
#======================================

# Precision
eps = 1e-6

#--------------------------------------
# Simple iteration method
#--------------------------------------

x_real = 0.5 # real value
x_start = 0  # start approximation

iterations = iterate(x_start, x_next_iteration_sim, eps)
build_graph(f"SIM: x(n+1) = 0,6 - 0,4 * (x(n))^2, x={x_real}", "SIM", iterations)

#--------------------------------------
# Newton's method
#--------------------------------------

x_real = -3  # real value
x_start = -2 # start approximation

iterations = iterate(x_start, x_next_iteration_newton, eps)
build_graph(f"SIM: x(n+1) = x(n) - f(x(n)) / f'(x(n)), x={x_real}", "NEWTON", iterations)
