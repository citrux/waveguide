# coding: utf-8
###############################################################################
# графическое решение дисперсионного уравнения
# совместимо и с python 3, и c python 3
###############################################################################

DEBUG = False

import matplotlib
import numpy as np
import matplotlib.pyplot as plt

from math import pi, tan

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

def bisection(f, left, right):
    """
    Метод бисекции поиска корня (трансцендентного) уравнения
    """
    precision = 0.005
    center = (left + right) / 2.0
    if right-left < precision:
        return center
    elif f(left) * f(right) > 0:
        return False
    elif abs(f(left)) < precision:
        return left
    elif abs(f(right)) < precision:
        return right
    else:
        if f(center) * f(left) <= 0:
            return bisection(f, left, center)
        else:
            return bisection(f, center, right)


def solution(condition, n1, n2):
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


def wavenumbers_relationship(omega):
    U1 = np.linspace(0, 7, 100)
    u2 = lambda u1: ((e2 * m2 - e1 * m1) * omega ** 2 / c ** 2 + u1 ** 2) ** 0.5
    U2 = list(map(u2, U1))
    plt.plot(U1, U2, "k--")


def plot(condition):
    for i in range(9):
        n = 2 * i + 1
        for i in range(n+1):
            u1, u2 = solution(condition, i, n - i)
            plt.plot(u1, u2, "k-")
            # нарисуем границы прямоугольников
            left = i * pi / 2.0 / l1
            right = (i + 1) * pi / 2.0 / l1
            bottom = (n - i) * pi / 2.0 / l2
            top = (n - i + 1) * pi / 2.0 / l2
            plt.plot([right, left, left, right, right],
                    [bottom, bottom, top, top, bottom], "k:")

    for omega in OMEGA:
        wavenumbers_relationship(omega)

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


if __name__ == '__main__':
    plot(e_condition)
    plot(m_condition)

