#include <cstdlib>
#include <iostream>

#include "numdiff.hpp"

long double calc_error(const Function_info& func_info, 
                        Num_der num_der, long double val, long double step) {
  return abs(num_der(func_info.func, val, step) - func_info.der(val));
}