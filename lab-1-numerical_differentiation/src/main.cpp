#include <cmath>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "numdiff.hpp"

using std::sin, std::cos, std::exp, std::sqrt, std::log;
using std::cerr, std::endl;

namespace {

// Numerical derivative formulas

const Num_der nd1 = 
[](Func f, long double x, long double h) -> long double { return (f(x + h) - f(x)) / h; };

const Num_der nd2 = 
[](Func f, long double x, long double h) -> long double { return (f(x) - f(x - h)) / h; };

const Num_der nd3 = 
[](Func f, long double x, long double h) -> long double { return (f(x + h) - f(x - h)) / (2.*h); };

const Num_der nd4 = 
[](Func f, long double x, long double h) -> long double {

  return (4. / 3.) * nd3(f, x,    h)  
       - (1. / 3.) * nd3(f, x, 2.*h);
};

const Num_der nd5 = 
[](Func f, long double x, long double h) -> long double {

  return (3. / 2.)  * nd3(f, x,    h) 
       - (3. / 5.)  * nd3(f, x, 2.*h) 
       + (1. / 10.) * nd3(f, x, 3.*h);
};

const std::array<Num_der, 5> Num_ders = {nd1, nd2, nd3, nd4, nd5};

const std::array<Function_info, 5> Comput_data = {{

  {
    .str_formula = "sin(x^2)",
    .func = [](long double x) -> long double { return sin(x * x); },
    .der  = [](long double x) -> long double { return cos(x * x) * 2. * x; }
  },

  {
    .str_formula = "cos(sin(x))",
    .func = [](long double x) -> long double { return cos(sin(x)); },
    .der  = [](long double x) -> long double { return - (sin(sin(x)) * cos(x)); }
  },

  {
    .str_formula = "exp(sin(cos))",
    .func = [](long double x) -> long double { return exp(sin(cos(x))); },
    .der  = [](long double x) -> long double { return - (exp(sin(cos(x))) 
                                   * cos(cos(x)) * sin(x)); }
  },

  {
    .str_formula = "ln(x+3)",
    .func = [](long double x) -> long double { return log(x + 3.); },
    .der  = [](long double x) -> long double { return 1. / (x + 3.); }
  },

  {
    .str_formula = "sqrt(x+3)",
    .func = [](long double x) -> long double { return sqrt(x + 3.); },
    .der  = [](long double x) -> long double { return 1. / (2. * sqrt(x + 3.)); }
  },
}};

}; // namespace

int main() {

  const std::string res_filename = "res/comput_results.txt";

  std::ofstream comput_results(res_filename, std::ios::out | std::ios::trunc);
  if (!comput_results.is_open()) {

    cerr << "Failed to open " << res_filename << endl;
    return EXIT_FAILURE;
  }

  const long double x_val = 2.28;
  comput_results << "x " << x_val << endl << endl;

  comput_results << std::setprecision(std::numeric_limits<double>::digits10);

  for (const auto& func_info : Comput_data) {

    comput_results << "name " << func_info.str_formula << endl;

    for (const auto& num_der : Num_ders) {

      comput_results << "data ";
      long double step = 1.;
      for (unsigned n = 1; n <= 21; n++) {

        comput_results << calc_error(func_info, num_der, x_val, step) << " ";
        step /= 2.;
      }

      comput_results << endl;
    }

    comput_results << endl;
  }

  comput_results.close();
  return EXIT_SUCCESS;
}