## Численное дифференцирования 

### Overview
В данной лабораторной работе предлагалось построить графики 
абсолютной погрешности различных методов численного дифферецирования в зависимости от шага для некоторых математических функций. 

### Функции
Список математических функций, на которых проводились вычисления:
$$1. y = sin(x^2)$$
$$2. y = cos(sin(x))$$
$$3. y = exp(sin(cos(x)))$$
$$4. y = ln(x+3)$$
$$5. y = sqrt(x+3)$$

## Формулы численного дифференцирования
Список предложенных формул для вычисления производной методом численного дифференцирования:
$$1. {f(x+h)-f(x) \over h}$$
$$2. {f(x)-f(x-h) \over h}$$
$$3. {f(x+h)-f(x-h) \over 2h}$$
$$4. {4 \over 3}{f(x+h)-f(x-h) \over 2h} - {1 \over 3}{f(x+2h)-f(x-2h) \over 4h}$$
$$5. {3 \over 2}{f(x+h)-f(x-h) \over 2h} - {3 \over 5}{f(x+2h)-f(x-2h) \over 4h} + {1 \over 10}{f(x+3h)-f(x-3h) \over 6h}$$

## Структура проекта
Проект включает в себя две части: программа на C++, предназначенная для выполнения необходимых вычислений, и <code>python</code>-скрипт для построения графиков с использованием <code>matplotlib pyplot</code>.

Заранее подготовленные графики находятся в директории <code>img/</code> в директории данной лабораторной работы. 

В следующем пункте предоставляется инструкция по сборке проекта и воспроизведению вычислений и построений графиков.

## Сборка и использование
Для того, чтобы собрать проект и воспроизвести результаты с построением графиков, необходимо:
 - Склонировать данный репозиторий и перейти в директорию <code>lab-1-numerical_differentiation</code>
   <code>git clone https://github.com/RustamSubkhankulov/computational-math.git</code>
   <code>cd lab-1-numerical_differentiation</code>
 - Воспользоваться <code>make plots</code>, чтобы собрать проект, провести вычисления и построить графики.

## Результаты
Далее представлены построенные программой графики. По обеим осям установлен логарифмический масштаб, каждый график содержит данные для всех формул численного дифференцирования для каждой функции.

Значение точки, в которой были произведены вычисления, выбрана произвольно.

<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-1-numerical_differentiation/img/0.png" alt="sin(x^2)" width="450"/>
<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-1-numerical_differentiation/img/1.png" alt="cos(sin(x))" width="450"/>
<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-1-numerical_differentiation/img/2.png" alt="exp(sin(cos(x)))" width="450"/>
<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-1-numerical_differentiation/img/3.png" alt="ln(x+3)" width="450"/>
<img src="https://github.com/RustamSubkhankulov/computational-math/blob/main/lab-1-numerical_differentiation/img/4.png" alt="exp(x+3)" width="450"/>

### Appendix. Прочие возможности системы сборки
Система сборки предоставляет ряд прочих средств:
 - <code>make build</code> - сборка проекта на С++ без запуска
 - <code>make run</code> - запуск вычислений с сборкой по необходимости без построения графиков
