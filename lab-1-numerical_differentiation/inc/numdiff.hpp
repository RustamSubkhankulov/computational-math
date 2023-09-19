#pragma once

#include <functional>
#include <utility>
#include <string>

// Function of one variable
using Func = std::function<long double(long double)>;

// Numerical derivative of function of one variable
using Num_der = std::function<long double(Func func, long double val, long double step)>;

struct Function_info {
  std::string str_formula; // string representation of formula
  Func func;               
  Func der;
};

long double calc_error(const Function_info& func_info, 
                        Num_der num_der, long double val, long double step);