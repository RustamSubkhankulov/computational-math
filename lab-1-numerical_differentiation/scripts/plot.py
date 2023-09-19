#!/usr/bin/python3

#======================================
# imports
#======================================

import sys
import matplotlib.pyplot  as plt

#======================================
# functions
#======================================

def read_comput_res(filename: str):
  with open(filename, "r") as comput_res:
    return [line.strip().split() for line in comput_res if line.strip()]

def get_graph_data(list: list):
  return [float(elem) for elem in list[1:]]

def calc_step(n: int) -> float:
  return 1. / 2 ** (n - 1)

def get_step_values(n: int):
  return [calc_step(i) for i in range(n)]

def build_graph(title: str, graph_data):

  plt.title(title)

  plt.xscale("log")
  plt.xlabel("ln(h)")

  plt.yscale("log")
  plt.ylabel("ln(error)")

  for num_der_idx in range(len(graph_data)):

    dots_num = len(graph_data[num_der_idx])
    step_values = get_step_values(dots_num)

    plt.plot(step_values, graph_data[num_der_idx], ".-")

  plt.savefig("img/graph_{}.png".format(title))
  plt.clf()

#======================================
# main 
#======================================

input = read_comput_res("res/comput_results.txt")

assert input[0][0] == "x"
x_value = float(input[0][1])

comput_res = input[1:]
line_idx = 0
line_num = len(comput_res)

while line_idx < line_num and comput_res[line_idx][0] == "name":

  str_formula = comput_res[line_idx][1]
  line_idx += 1

  graph_data = []

  while line_idx < line_num and comput_res[line_idx][0] == "data":

    graph_data.append(get_graph_data(comput_res[line_idx]))
    line_idx += 1

  build_graph(str_formula + " at x = {}".format(x_value), graph_data)
