## Интерполяция

### Overview
В данной лабораторной работе предлагалось интерполировать зависимость населения США от года по десяти исходным точкам - от 1910 до 2000 и сравнить полученное по экстраполяции значение населения в 2010 году с реальным. 

```
h = 10.
x = [1910. + h * i for i in range(10)]
y = [ 92228496., 106021537., 123202624., 132164569., 151325798., 
     179323175., 203211926., 226545805., 248709873., 281421906. ]

amount = len(y)
data = {x[i] : y[i] for i in range(amount)}

x0 = 2010
y_real = 308745538
```

### Методы решения поставленной задачи:
1. Интерполянт в форме Ньютона  
2. Сплайн экстраполяция
3. Метод Наименьших Квадратов

### Структура проекта
Проект включает в себя <code>python</code>-скрипт для проведения вычислений и построения графиков с использованием <code>matplotlib pyplot</code>.

Заранее подготовленные графики и проие изображения, включенные в отчёт, находятся в директории <code>img/</code> в директории данной лабораторной работы. 

В следующем пункте предоставляется инструкция по запуску проекта и воспроизведению вычислений и построений графиков.

### Использование
Для того, чтобы воспроизвести результаты с построением графиков, необходимо:
 - Склонировать данный репозиторий и перейти в директорию <code>lab-2-SoLAE</code>
   1. <code>git clone https://github.com/RustamSubkhankulov/computational-math.git</code>
   2. <code>cd lab-4-interpolation</code>
 - Воспользоваться <code>./scripts/interpolation.py</code> для интерполяции зависимости населения от года, чтобы провести вычисления и построить графики. 

Все выходные файлы сохраняются в директории <code>res/</code> относительно директории работы. 

### Результаты

Ниже приведены графики полученных интерполянтов для различных методов

1. Интерполянт в форме Ньютона:

```
def div_diff(x, y, k: int, j: int) -> float:
  """
  Calculate divided difference in x_{j} of order k.
  """

  if k == 0:
    return y[j]

  return (div_diff(x, y, k-1, j+1) - div_diff(x, y, k-1, j)) / (x[j + k] - x[j])
```

```
def build_interpolant_newton(x, y):
  """
  Calculate coeffs of interpolant in Newtonian form.
  """

  coeffs = []
  for k in range(len(x)):
    coeffs.append(div_diff(x, y, k, 0))

  return coeffs

def calculate_newton(x, coeffs, x0) -> float:
  """
  Calculate value of interpolant in Newtonian form in x0.
  """

  value = 0.
  for i in range(len(coeffs)):

    term = coeffs[i]

    for j in range(i):
      term *= (x0 - x[j])

    value += term

  return value
```

<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-4-interpolation/img/newton.png" alt="newton" width="450"/>

Результаты экстраполяции:
```
Interpolant in Newtonian form:
In point 2010: Res == 827906509; Real == 308745538;
Diff == 519160971;
```

Может видеть, что полученное значение сильно отличатеся от реального. Это означает, что метод Ньютона плохо подходит для данной задачи.

2. Сплайн экстраполяция:

```
def diff(x, y, first, last):
  """
  Helper for calculating right side.
  """ 

  if first == last:
      return y[first]
  return (diff(x, y, first+1, last) - diff(x, y, first, last-1)) / (x[last] - x[first])

def u(x, y, i):
  """
  Helper for calculating right side.
  """  

  return (diff(x, y, i, i + 1) - diff(x, y, i - 1, i)) / (x[i + 1] - x[i - 1])

def f_elem_spline(x, y, order: int, i: int) -> float:
  """
  Returs ith element of right side of the system.
  """

  return 6. * u(x, y, i)

def get_f_spline(x, y, order: int):
  """
  Returs right side of the system.
  """

  return [f_elem_spline(x, y, order, i) for i in range(order)]

def matrix_elem_spline(h: float, x, order: int, i: int, j: int) -> float:
  """
  Returns (i,j) matrix element.
  """

  if i == j:
    return 2.

  if abs(i - j) > 1:
    return 0
  else:
    return 0.5

def get_matrix_spline(h: float, x, y, order: int):
  """
  Returns system's matrix.
  """

  matrix = np.empty((order, order))
  for i in range(order):
    for j in range(order):
      matrix[i][j] = matrix_elem_spline(h, x, order, i, j)

  return matrix

def build_interpolant_spline(h: float, x, y):
  """
  Calculate coeffs of spline approximation.
  """

  order = len(x) - 2

  matrix = get_matrix_spline(h, x, y, order)
  f = get_f_spline(x, y, order)

  print([0, *gauss(matrix, f), 0])
  return [0, *gauss(matrix, f), 0]

def get_spline_index(x, x0: float) -> int:
  """
  Get index of spline for target point.
  """

  ind = 0

  for i in range(len(x)):
    if x0 <= x[i]:
      ind = i
      break

  return ind

def calculate_spline(h: float, x, y, coeffs, x0) -> float:
  """
  Calculate value of spline approximation in x0.
  """
  
  ind = get_spline_index(x, x0)
  x_i = x[ind]

  a = y[ind]
  b = coeffs[ind] * h / 3 + coeffs[ind-1] * h / 6 + diff(x, y, ind-1, ind)
  c = coeffs[ind]
  d = (coeffs[ind] - coeffs[ind-1]) / h

  return a + b * (x0-x_i) + (c/2) * (x0-x_i)**2 + (d/6) * (x0-x_i)**3
```

<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-4-interpolation/img/spline.png" alt="spline" width="450"/>

```
Spline approximation:
In point 2010: Res == 302443396.0; Real == 308745538;
Diff == 6302142;
```

Полученное значение разности мало, и это показывает, что зависимость хорошо экстраполируется методов кубической сплайн экстраполяции.

3. Метод Наименьших Квадратов:

```
def f_elem_lsm(x, y, order: int, i: int) -> float:
  """
  Returs ith element of right side of the system.
  """

  power = order - (i + 1)
  return sum([y[k] * pow(x[k], power) for k in range(len(x))])

def get_f_lsm(x, y, order: int):
  """
  Returs right side of the system.
  """

  return [f_elem_lsm(x, y, order, i) for i in range(order)]

def matrix_elem_lsm(x, order: int, i: int, j: int) -> float:
  """
  Returns (i,j) matrix element.
  """

  power = 2 * order - (i + j + 2)
  return sum(pow(x[k], power) for k in range(len(x)))
  
def get_matrix_lsm(x, order: int):
  """
  Returns system's matrix.
  """

  matrix = np.empty((order, order))
  for i in range(order):
    for j in range(order):
      matrix[i][j] = matrix_elem_lsm(x, order, i, j)

  return matrix

def build_interpolant_lsm(x, y, order: int):
  """
  Calculate coeffs of least square method interpolation.
  """

  matrix = get_matrix_lsm(x, order)
  f = get_f_lsm(x, y, order)

  return gauss(matrix, f)

def calculate_lsm(coeffs, x0) -> float:
  """
  Calculate value of least square method interpolant in x0.
  """

  return sum([coeffs[i] * pow(x0, len(coeffs) - (i + 1)) for i in range(len(coeffs))])
```

<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-4-interpolation/img/lsm.png" alt="lsm" width="450"/>

```
Least square method:
In point 2010: Res == 312470336.42485046; Real == 308745538;
Diff == 3724798;
```
Зависимость населения очень хорошо получилось приблизить квадратичным полиномом с помощью Метода Наименьших квадратов.
