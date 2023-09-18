#include <cstdlib>

#include "numdiff.hpp"

double calc_error(const Function_info& func_info, 
                        Num_der num_der, double val, double step) {
  return abs(num_der(func_info.func, val, step) - func_info.der(val));
}