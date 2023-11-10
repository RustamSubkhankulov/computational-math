## Решение НАУ и СНАУ

### Overview
В данной лабораторной работе предлагалось решить нелинейное алгебраическое уравнения и систему нелинейных алгебраических уравнений. 
Чтобы показать, что найденное решение является верным построены графики сходимости к корням.

НАУ:
$$2x^2 + 5x - 3 = 0$$

CНАУ:
$$sin(x + 1) - y = 1.5$$
$$x + cos(y - 2) = 0.5$$

### Методы решения НАУ:
1. Метод простой итерации для нелинейного алгебраического уравнения:
  $$x(n+1) = 0,6 - 0,4 * (x(n))^2$$  
2. Метод Ньютона для нелинейного алгебраического уравнения:
  $$x(n+1) = x(n) - f(x(n)) / f'(x(n))$$ 

### Методы решения СНАУ:
1. Метод простой итерации для систем нелинейных алгебраических уравнений
  $$x(n+1) = 0.5 - cos(y(n) - 2)$$
  $$y(n+1) = sin(x(n) + 2) - 1.5$$
2. Метод Ньютона для систем нелинейных алгебраических уравнений 
  $$\ves{x(n+1)} = \ves{x(n)} - J^(-1)(\ves{x(n)}) * \ves{F(\ves{x(n)})}$$

### Структура проекта
Проект включает в себя <code>python</code>-скрипт для проведения вычислений и построения графиков с использованием <code>matplotlib pyplot</code>.

Заранее подготовленные графики и проие изображения, включенные в отчёт, находятся в директории <code>img/</code> в директории данной лабораторной работы. 

В следующем пункте предоставляется инструкция по запуску проекта и воспроизведению вычислений и построений графиков.

### Использование
Для того, чтобы воспроизвести результаты с построением графиков, необходимо:
 - Склонировать данный репозиторий и перейти в директорию <code>lab-2-SoLAE</code>
   1. <code>git clone https://github.com/RustamSubkhankulov/computational-math.git</code>
   2. <code>cd lab-3-SoNLAE</code>
 - Воспользоваться <code>./scripts/nlae.py</code> для решения нелинейного алгебраического уравнения или <code>./scripts/sonlae.py</code> для решения системы нелинейных алгебраических уравнений, чтобы провести вычисления и построить графики. 

Все выходные файлы сохраняются в директории <code>res/</code> относительно директории работы. 

### Результаты

Ниже приведены графики сходимости последовательных итераций к корням уравнений и систем уравнений.

#### НАУ

Общий метод итерации:
```
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
```

Метод простой итерации:
```
def x_next_iteration_sim(x_n: float) -> float:
  """
  x(n+1) = 0,6 - 0,4 * (x(n))^2
  """
  return 0.6 - 0.4 * pow(x_n, 2.0)
```

<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-3-SoNAE/img/SIM.png" alt="SIM" width="450"/>

Метод Ньютона:
```
def x_next_iteration_newton(x_n: float) -> float:
  """
  x(n+1) = x(n) - f(x(n)) / f'(x(n))
  """
  return x_n - (f(x_n)) / (f_deriv(x_n))
```

<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-3-SoNAE/img/NEWTON.png" alt="NEWTON" width="450"/>

#### СНАУ

Общий метод итерации:
```
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
```

Метод простой итерации:
```
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
```

<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-3-SoNAE/img/SIM_SYS.png" alt="SIM_SYS" width="450"/>

Метод Ньютона:
```
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
```

<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-3-SoNAE/img/NEWTON_SYS.png" alt="NEWTON_SYS" width="450"/>
