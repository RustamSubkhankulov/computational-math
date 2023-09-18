#pragma once

#include <functional>
#include <utility>
#include <string>

// Function of one variable
using Func = std::function<double(double)>;

// Numerical derivative of function of one variable
using Num_der = std::function<double(Func func, double val, double step)>;

struct Function_info {
  std::string str_formula; // string representation of formula
  Func func;               
  Func der;
};

double calc_error(const Function_info& func_info, 
                        Num_der num_der, double val, double step);