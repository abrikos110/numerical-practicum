import os
import subprocess as sub

import numpy
import matplotlib.pyplot as plt

from numpy import exp, sqrt

c = sqrt(2/5)
cc = 4 - sqrt(10) + (4 + sqrt(10)) * exp(2*c)
g = lambda x: exp(-0.5) * 2 * (cc - 4 * (exp(c * (1+x)) + exp(c * (1-x)))) / cc

os.system('g++ main.cpp -o main -Wall -Wextra -Wpedantic -O3 -std=c++17')
s = sub.getoutput('echo ' + input() + ' | ./main').split('\n')

m = [numpy.array(eval('[' + s[i] + ']')) for i in range(-5, -2)]
v = numpy.array(eval('[' + s[-1] + ']'))
f = numpy.array(eval('[' + s[-2] + ']'))

xx = numpy.linspace(0, 1, len(v))
y = g(xx)
# m @ y - f

my = y * m[1]
my[:-1] += m[0] * y[1:]
my[1:] += m[2] * y[:-1]

print('approx:', max(abs(my)))
print('error:', max(abs(y - v)))

plt.plot(xx, v, label='calc')
plt.plot(xx, y, label='orig')
plt.legend()
plt.show()
