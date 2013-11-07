# coding: utf-8

import matplotlib
import numpy as np
import matplotlib.pyplot as plt

from math import pi, tan

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

e1 = 1.0
e2 = 5.0
l1 = 3.5
l2 = 1.5

OMEGA = np.linspace(2e10,6e10,5)
c = 3e10

def fun(k1, k2):
    return k1 * tan(k1 * l1) / e1 + k2 * tan(k2 * l2) / e2

def dichotomy(f, a, b):
    eps = 0.005
    c = (a+b)/2
    if b-a < eps:
        return c
    elif f(a)*f(b) > 0:
        return False
    elif abs(f(a)) < eps:
        return a
    elif abs(f(b)) < eps:
        return b
    else:
        if f(c)*f(a) <= 0:
            return dichotomy(f, a, c)
        else:
            return dichotomy(f, c, b)


def solution(n1, n2):
    delta = 0.005
    K1 = np.arange(n1*pi/2/l1+delta, (n1+1)*pi/2/l1-delta, delta)
    x,y = [],[]
    for i in K1:
        j = dichotomy(lambda x: fun(i,x),
                n2*pi/2/l2+delta, (n2+1)*pi/2/l2-delta)
        if j:
            x.append(i)
            y.append(j)
    plt.plot(x,y,"k-")

for i in range(9):
    n = 2 * i + 1
    for i in range(n+1):
        solution(i,n-i)

def hyperbole(omega):
    k1_list = np.linspace(0, 4, 100)
    k2_list = list(map(lambda x: ((e2-e1)*omega**2/c**2 + x**2), k1_list))
    plt.plot(k1_list, k2_list, "k--")

for omega in OMEGA:
    hyperbole(omega)

plt.xlabel(r"$k_1, cm^{-1}$")
plt.ylabel(r"$k_2, cm^{-1}$")

plt.show()
