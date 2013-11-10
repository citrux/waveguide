# coding: utf-8
###############################################################################
# графическое решение дисперсионного уравнения
# совместимо и с python 3, и c python 3
###############################################################################

DEBUG = False

import matplotlib
import numpy as np
import matplotlib.pyplot as plt

from math import pi, tan, log

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

e1, m1, l1 = 1.0, 1.0, 3.5
e2, m2, l2 = 5.0, 1.0, 1.5
# скорость света, см/с
c = 3e10
# частоты (СВЧ 3 ÷ 30 ГГц) → omega = 20 ÷ 200 Град/с
OMEGA = np.linspace(2e10, 6e10, 5)

def e_condition(u1, u2):
    return u1 * tan(u1 * l1) / e1 + u2 * tan(u2 * l2) / e2

def m_condition(u1, u2):
    return m1 * tan(u1 * l1) / u1 + m2 * tan(u2 * l2) / u2

def bisection(f, left, right, precision=5e-3):
    """
    Метод бисекции поиска корня (трансцендентного) уравнения
    """
    center = (left + right) / 2.0
    if right-left < precision:
        return center
    elif f(left) * f(right) > 0:
        return False
    #elif abs(f(left)) < precision:
        #return left
    #elif abs(f(right)) < precision:
        #return right
    else:
        if f(center) * f(left) <= 0:
            return bisection(f, left, center)
        else:
            return bisection(f, center, right)


def curve(condition, n1, n2):
    delta = 0.005

    # границы области, в которой ищется решение
    left = n1 * pi / 2.0 / l1
    right = (n1 + 1) * pi / 2.0 / l1
    bottom = n2 * pi / 2.0 / l2
    top = (n2 + 1) * pi/2.0/l2

    U1 = np.arange(left + delta, right - delta, delta)
    x,y = [],[]
    for u1 in U1:
        u2 = bisection(lambda u2: condition(u1, u2), bottom + delta, top - delta)
        if u2:
            x.append(u1)
            y.append(u2)
    return x, y


def plot_wavenumbers_relationship(omega, color="k", right=7):
    U1 = np.linspace(0, right, 100)
    u2 = lambda u1: ((e2 * m2 - e1 * m1) * omega ** 2 / c ** 2 + u1 ** 2) ** 0.5
    U2 = list(map(u2, U1))
    plt.plot(U1, U2, "%s-" % color)


def solution(condition, n1, n2, omega, delta=1e-4):
    test = lambda u1, u2: (u1 ** 2 - u2 ** 2 -\
                    omega ** 2 / c ** 2 * (e1 * m1 - e2 * m2))
    # u1
    left = n1 * pi / 2.0 / l1 + delta
    right = (n1 + 1) * pi / 2.0 / l1 - delta
    # u2
    bottom = n2 * pi / 2.0 / l2 + delta ** 2    # немного магии: из-за "плохих"
    top = (n2 + 1) * pi / 2.0 / l2 - delta ** 2 # производных приходится ставить
                                                # разные границы
    u2 = lambda u1: bisection(lambda u2: condition(u1, u2), bottom, top, delta)
    result = False
    sleft = u2(left)
    sright = u2(right)

    #print(test(left, sleft) * test(right, sright))
    while test(left, sleft) * test(right, sright) <= 0:
        #print(left, right, sleft, sright)
        center = (left + right) / 2.0
        if test(left, sleft) * test(center, u2(center)) <= 0:
            right = center
            sright = u2(center)
        elif (right-left) <= delta:
            if abs(test(left, sleft)) < abs(test(right, sright)):
                result = (left, sleft)
                break
            else:
                result = (right, sright)
                break
        else:
            left = center
            sleft = u2(center)
    return result

def transversal_wavenumbers(condition, n, omega, delta):
    for i in range(2*n):
        a = solution(condition, i, 2*n - 1 - i, omega, delta)
        if a:
            u1, u2 = a
            b = -log(delta) / log(10)
            u1 = round(u1 * 10 ** b) * 10 ** -b
            u2 = round(u2 * 10 ** b) * 10 ** -b
            return u1, u2
    return False, False

def plot_transversal(condition):
    m = 3
    k = 0
    color = ["r", "y", "g", "c", "b", "m"]
    for omega in OMEGA:
        plot_wavenumbers_relationship(omega, color[k], (m + 1) * pi / 2.0 / l1)
        k += 1
    for i in range(m):
        n = 2 * i + 1
        for j in range(n+1):
            U1, U2 = curve(condition, j, n - j)
            # нарисуем границы прямоугольников
            left = j * pi / 2.0 / l1
            right = (j + 1) * pi / 2.0 / l1
            bottom = (n - j) * pi / 2.0 / l2
            top = (n - j + 1) * pi / 2.0 / l2
            plt.plot([right, left, left, right, right],
                    [bottom, bottom, top, top, bottom], "k:")
            plt.plot(U1, U2, "k-")

        for omega in OMEGA:
            u1, u2 = transversal_wavenumbers(condition, i + 1, omega, 1e-7)
            plt.plot([u1], [u2], "ro")

    plt.xlabel(r"$u_1, cm^{-1}$")
    plt.ylabel(r"$u_2, cm^{-1}$")

    if DEBUG:
        plt.show()
    else:
        if condition is e_condition:
            name = "e"
        elif condition is m_condition:
            name = "m"
        else:
            name = "wtf"
        plt.savefig(name + ".png")
    plt.cla()


def longitudinal_wavenumber(condition, n):
    H, O = [], []
    for omega in np.linspace(2e10, 6e10, 100):
        u1, u2 = transversal_wavenumbers(condition, n, omega, 1e-7)
        if u1 and u2 and (omega/c * (e1 * m1) ** 0.5 > u1):
            O.append(omega)
            H.append((omega ** 2 / c**2 * e1 * m1 - u1 ** 2) ** 0.5)
    return O, H

def plot_longitudinal(condition, N):
    for i in N:
        O, H = longitudinal_wavenumber(condition, i)
        plt.plot(O, H, "k-")

    plt.xlabel(r"$\omega, rad/s$")
    plt.ylabel(r"$h, cm^{-1}$")
    if DEBUG:
        plt.show()
    else:
        if condition is e_condition:
            name = "e_h"
        elif condition is m_condition:
            name = "m_h"
        else:
            name = "wtf"
        plt.savefig(name + ".png")
    plt.cla()

if __name__ == '__main__':
    plot_longitudinal(e_condition, [1,2,3])
    plot_longitudinal(m_condition, [2,3])




