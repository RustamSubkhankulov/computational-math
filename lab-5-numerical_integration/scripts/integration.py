#!/usr/bin/python3

# Calculation of the definite integral
# by the trapezoidal method with 
# and without Richardson extrapolation 
# and Simpson's method

#======================================
# Imports
#======================================

import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from os import mkdir

#======================================
# Methods
#======================================

#--------------------------------------
# Trapezoidal method
#--------------------------------------

def num_integral_trapezoidal_cotes(h: float, y: list[float]) -> float:
  """
  Calculate integral using Trapezoidal method via Cotes formula
  """
  return h * ((y[0] + y[len(y) - 1]) / 2 + sum(y[1:len(y) - 1]))

#--------------------------------------
# Simpson's method
#--------------------------------------

def num_integral_simpsons(h: float, x: list[float], y: list[float]) -> float:
  """
  Calculate integral using Simpson's method
  """

  n = len(y) - 1
  return (y[0] + y[n] + sum([(2 + 2 * (i % 2 == 1)) * y[i] for i in range(1, n)])) * h / 3

#--------------------------------------
# Misc
#--------------------------------------

def print_res(title: str, answer: float, got: float):
  
  print(title)
  print(f"Calculated result: {got:.6f}; Answer: {answer:.6f};")

#======================================
# Main 
#======================================

def main():

  #---------------
  # Table function
  #---------------
  amount = 9
  h = 0.25

  x = [i * h for i in range(amount)]

  y = [ 1.000000, 0.989616, 0.958851,
        0.908852, 0.841471, 0.759188,
        0.664997, 0.562278, 0.454649 ]

  data = {x[i] : y[i] for i in range(amount)}

  # Real calculated answer 
  I_answer = 1.605412976802695

  #--------------------------
  # Trapezoidal method h=0.25
  #--------------------------
  Ih = num_integral_trapezoidal_cotes(h, y)
  print_res("Trapezoidal method h=0.25", I_answer, Ih)

  #-------------------------
  # Trapezoidal method h=0.5
  #-------------------------
  I2h = num_integral_trapezoidal_cotes(2*h, [y[i] for i in range(0, len(y), 2)])
  print_res("Trapezoidal method h=0.5", I_answer, I2h)

  #-------------------------
  # Richardson extrapolation
  #-------------------------
  Ir = Ih + (Ih - I2h) / (pow(2,2) - 1)
  print_res("Richardson extrapolation", I_answer, Ir)

  #-----------------
  # Simpson's method
  #-----------------
  Is = num_integral_simpsons(h, x, y)
  print_res("Simpson's method", I_answer, Is)

#======================================

if __name__ == '__main__':
  main()
