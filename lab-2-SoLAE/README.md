## Решение СЛАУ прямыми и итерационными методами 

### Overview
В данной лабораторной работе предлагалось решить СЛАУ некоторыми прямыми и итерационными методами. 
Чтобы показать, что найденное решение является верным 
- для прямых методов: подставить решение в исходную систему и получить невязку ~ 0. 
- для итерационных методов: показать убывание невязки (график изменения невязки от числа итераций).

<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-2-SoLAE/img/system.png" alt="system" width="450">

### Прямые методы решения СЛАУ:
1. Гаусса с выбором главного элемента
2. LU-разложение

### Итерационные методы решения СЛАУ:
1. Метод Зейделя
2. Метод Якоби
3. Метод верхней релаксации

### Структура проекта
Проект включает в себя <code>python</code>-скрипт для проведения вычислений и построения графиков с использованием <code>matplotlib pyplot</code>.

Заранее подготовленные графики и проие изображения, включенные в отчёт, находятся в директории <code>img/</code> в директории данной лабораторной работы. 

В следующем пункте предоставляется инструкция по запуску проекта и воспроизведению вычислений и построений графиков.

### Использование
Для того, чтобы воспроизвести результаты с построением графиков, необходимо:
 - Склонировать данный репозиторий и перейти в директорию <code>lab-2-SoLAE</code>
   1. <code>git clone https://github.com/RustamSubkhankulov/computational-math.git</code>
   2. <code>cd lab-2-SoLAE</code>
 - Воспользоваться <code>./scripts/solae.py</code>, чтобы провести вычисления и построить графики. 

Все выходные файлы сохраняются в директории <code>res/</code> относительно директории работы. 

### Результаты

Для определения совпадения результатов решения системы прямыми и итерационными методами было взято значение точности, равное 10^-6
```{python}
eps = 1e-6
```
Норма вектора:
```{python}
def norm(vec):
  return np.max(np.abs(vec))
```

#### Прямые методы

1. Метод Гаусса.
Непосредственно сам метод:
```{python}
def gauss(matrix, f):

  ext = get_extended_matrix(matrix, f)  

  n = len(matrix)
  for i in range(n):

      leading = i + np.argmax(np.abs(matrix[:,i][i:]))
      ext[[i, leading]] = ext[[leading, i]] 

      ext[i] /= ext[i][i]
      row = ext[i]

      for r in ext[i + 1:]:
          r -= r[i] * row

  for i in range(n - 1, 0, -1):
      row = ext[i]
      for r in reversed(ext[:i]):
          r -= r[i] * row

  return ext[:,-1]
```
Вычисление решения и проверка результатов:
```{python}
x = gauss(matrix, f)
assert(norm(np.matmul(matrix,x) - f) < eps)
```
2. LU-разложение
Разложение матрицы:
```{python}
def lu_decomposition(matrix):
  # permutation, lower, upper = sp.linalg.lu(matrix)
  # return lower, upper

  n = len(matrix)

  lower = np.full((n, n), 0.)
  upper = np.full((n, n), 0.)

  for i in range(n):
    lower[i][i] = 1

    if i != 0:
      for j in range(i, n):
        upper[i][j] = matrix[i][j]
        for k in range(0, i):
            upper[i][j] -= lower[i][k] * upper[k][j]

      for j in range(i + 1, n):
        lower[j][i] = matrix[j][i]
        for k in range(0, i):
          lower[j][i] -= lower[j][k] * upper[k][i]
        lower[j][i] /= upper[i][i]
    else:
        for j in range(0, n):
            upper[0, j] = matrix[0, j]
            lower[j, 0] = matrix[j, 0] / matrix[0, 0]

  return lower, upper
```
Решение системы соответствующим методом:
```{python}
def lu(matrix, f):

  n = len(matrix)
  lower, upper = lu_decomposition(matrix)

  x = [0 for i in range(n)]
  y = [0 for i in range(n)]

  for i in range(1, n + 1):

    summ = 0
    for k in range(1, i):
      summ += lower[i-1][k-1] * y[k-1]

    y[i-1] = f[i-1] - summ

  for i in range(n, 0, -1):

    summ = 0
    for k in range(n, i, -1):
      summ += upper[i-1][k-1] * x[k-1]

    x[i-1] = (y[i-1] - summ) / upper[i-1][i-1]

  return x;
```
Проверка:
```{python}
x = lu(matrix, f)
assert(norm(np.matmul(matrix,x) - f) < eps)
```
#### Итерационные методы

Разложение на lower, diagonal и upper матрицы:
```{python}
def ldu_decomposition(matrix):
  u = np.triu(matrix)
  l = np.tril(matrix)

  lower = matrix - u
  diagonal = l + u - matrix
  upper = matrix - l

  return lower, diagonal, upper
```

Метод простой итерации:
```{python}
def simple_iteration_method(B, F, matrix, f, iterations_number):

  n = len(F)

  x = [0 for i in range(n)]
  err = [0 for i in range(iterations_number)]

  for iter in range(iterations_number):
    
    x = np.matmul(B, x) + F
    
    diff =  f - np.matmul(matrix, x)
    err[iter] = norm(diff) 

  return x, err)
```

1. Метод Зейделя
```{python}
def seidel(matrix, f, iterations_number):

  lower, diagonal, upper = ldu_decomposition(matrix)

  B = - np.matmul(np.linalg.inv(lower + diagonal), upper)
  F =   np.matmul(np.linalg.inv(lower + diagonal), f)

  return simple_iteration_method(B, F, matrix, f, iterations_number)
```
2. Метод Якоби
```{python}
def jacobi(matrix, f, iterations_number):

  lower, diagonal, upper = ldu_decomposition(matrix)

  B = - np.matmul(np.linalg.inv(diagonal), lower + upper)
  F =   np.matmul(np.linalg.inv(diagonal), f)

  return simple_iteration_method(B, F, matrix, f, iterations_number)
```
3. Метод верхней релаксации
```{python}
def successive_over_relaxation(matrix, f, iterations_number, w):

  lower, diagonal, upper = ldu_decomposition(matrix)

  B = - np.matmul(np.linalg.inv(diagonal + w * lower), (w - 1) * diagonal + w * upper)
  F =   np.matmul(np.linalg.inv(diagonal + w * lower), f) * w

  return simple_iteration_method(B, F, matrix, f, iterations_number)
```
