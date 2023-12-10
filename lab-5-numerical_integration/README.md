## Численное интегрирование

### Overview
В данной лабораторной работе предлагалось для функции, заданной таблично, вычислить значение определенного интеграла методом трапеций, сделать уточнение результата экстраполяцией Ричардсона. Сравнить уточненный результат с вычислениями по методу Симпсона.

| x_{i} | f(x_{i}) |
|-------|----------|
| 0.00  | 1.000000 |
| 0.25  | 0.989616 |
| 0.50  | 0.958851 |
| 0.75  | 0.908852 |
| 1.00  | 0.841471 |
| 1.25  | 0.759188 |
| 1.50  | 0.664997 |
| 1.75  | 0.562278 |
| 2.00  | 0.454649 |

```
amount = 9
h = 0.25

x = [i * h for i in range(amount)]

y = [ 1.000000, 0.989616, 0.958851,
      0.908852, 0.841471, 0.759188,
      0.664997, 0.562278, 0.454649 ]

data = {x[i] : y[i] for i in range(amount)}

# Real calculated answer 
I_answer = 1.605412976802695
```

### Методы численного интегрирования:
1. Метод трапеций с уточнением экстраполяцией Ричардсона 
2. Метод Симпсона

Во всех случаях применяется равномерная сетка.

### Структура проекта
Проект включает в себя <code>python</code>-скрипт для проведения вычислений.

В следующем пункте предоставляется инструкция по запуску проекта и воспроизведению вычислений.

### Использование
Для того, чтобы воспроизвести результаты, необходимо:
 - Склонировать данный репозиторий и перейти в директорию <code>lab-5-integration</code>
   1. <code>git clone https://github.com/RustamSubkhankulov/computational-math.git</code>
   2. <code>cd lab-5-integration</code>
 - Воспользоваться <code>./scripts/integration.py</code> для применения методов численного интегрирования для вычисления результатов. 

Пример вывода программы:
```
Trapezoidal method h=0.25
Calculated result: 1.603144; Answer: 1.605413;
Trapezoidal method h=0.5
Calculated result: 1.596322; Answer: 1.605413;
Richardson extrapolation
Calculated result: 1.605419; Answer: 1.605413;
Simpson's method
Calculated result: 1.605419; Answer: 1.605413;
```

### Результаты

1. Метод трапеций:

```
def num_integral_trapezoidal_cotes(h: float, y: list[float]) -> float:
  """
  Calculate integral using Trapezoidal method via Cotes formula
  """
  return h * ((y[0] + y[len(y) - 1]) / 2 + sum(y[1:len(y) - 1]))
```

Результаты расчетов для h=0,25 и h=0,5:
```
Trapezoidal method h=0.25
Calculated result: 1.603144; Answer: 1.605413;
Trapezoidal method h=0.5
Calculated result: 1.596322; Answer: 1.605413;
```

Вычисление интеграла с меньшим значением шага закономерно даёт значение ближе к реальному ответу. 

Эстраполяция Ричардсона:
```
Ir = Ih + (Ih - I2h) / (pow(2,2) - 1)
```

Результат вычислений:
```
Richardson extrapolation
Calculated result: 1.605419; Answer: 1.605413;
```

Применение эстраполяции Ричардсона приближает полученное из вычислений значение ещё ближе к ответу.

2. Метод Симпсона:

```
def num_integral_simpsons(h: float, x: list[float], y: list[float]) -> float:
  """
  Calculate integral using Simpson's method
  """

  n = len(y) - 1
  return (y[0] + y[n] + sum([(2 + 2 * (i % 2 == 1)) * y[i] for i in range(1, n)])) * h / 3
```

Результат вычислений:
```
Simpson's method
Calculated result: 1.605419; Answer: 1.605413;
```

Метод Симпсона показывает меньшее значение отклонения от реального ответа, чем метод трапеций.
