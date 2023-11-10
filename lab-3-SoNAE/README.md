## Решение НАУ и СНАУ

### Overview
В данной лабораторной работе предлагалось решить нелинейное алгебраическое уравнения и систему нелинейных алгебраических уравнений. 
Чтобы показать, что найденное решение является верным построены графики сходимости к корням.

НАУ:
$$2x^2 + 5x - 3 = 0$$

CНАУ:


### Методы решения НАУ:
1. Метод простой итерации для нелинейного алгебраического уравнения:
  $$x(n+1) = 0,6 - 0,4 * (x(n))^2$$  
2. Метод Ньютона для нелинейного алгебраического уравнения:
  $$x(n+1) = x(n) - f(x(n)) / f'(x(n))$$ 

### Методы решения СНАУ:
1. Метод простой итерации для систем нелинейных алгебраических уравнений
2. Метод Ньютона для систем нелинейных алгебраических уравнений 

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

<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-1-numerical_differentiation/img/NEWTON.png" alt="NEWTON" width="450"/>

#### СНАУ


