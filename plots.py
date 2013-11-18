# coding: utf-8
###############################################################################
# графическое решение дисперсионного уравнения
# совместимо и с python 3, и c python 3
# частоты (СВЧ 3 ÷ 30 ГГц) → omega = 20 ÷ 200 Град/с
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
b = 3.0
# скорость света, см/с
c = 3e10

def e_relation(u1, u2):
    return u1 * tan(u1 * l1) / e1 + u2 * tan(u2 * l2) / e2

def m_relation(u1, u2):
    return m1 * tan(u1 * l1) / u1 + m2 * tan(u2 * l2) / u2

def bisection(f, left, right, precision=5e-3):
    """
    Метод бисекции поиска корня (трансцендентного) уравнения
    """
    def sgn(x):
        if x > 0:
            return 1
        if x < 0:
            return -1
        return 0

    center = (left + right) / 2.0
    if right-left < precision:
        return center
    if sgn(f(left)) * sgn(f(right)) > 0:
        return False
    if sgn(f(center)) * sgn(f(left)) <= 0:
        return bisection(f, left, center)
    return bisection(f, center, right)


def transversal_wavenumbers(relation, m, omega, precision):
    for i in range(2*m):
        lt = ((2*m - i) / l2) ** 2 - (i / l1) ** 2 -\
                (2 * omega / c / pi) ** 2 * (e2 * m2 - e1 * m1)
        rd = ((2*m - i - 1) / l2) ** 2 - ((i + 1) / l1) ** 2 -\
                (2 * omega / c / pi) ** 2 * (e2 * m2 - e1 * m1)

        if lt * rd < 0:
            first = lambda x: x
            second = lambda x:\
                    (x ** 2 + (omega / c) ** 2 * (e2 * m2 - e1 * m1)) ** 0.5

            left = i * pi / 2.0 / l1
            right = (i + 1) * pi / 2.0 / l1
            down = (2 * m - i - 1) * pi / 2.0 / l2
            up = (2 * m - i) * pi / 2.0 / l2

            new_left = down ** 2 -  (omega / c) ** 2 * (e2 * m2 - e1 * m1)
            new_right = up ** 2 -  (omega / c) ** 2 * (e2 * m2 - e1 * m1)
            if new_left > left ** 2:
                left = new_left ** 0.5
            if new_right < right ** 2:
                right = new_right ** 0.5

            margin = 1e-7
            if relation(left + margin, second(left + margin)) <= 0:
                u1 = bisection(lambda x: relation(first(x),\
                    second(x)), left + margin, right - margin, precision)
                u2 = second(u1)

                if DEBUG:
                    print("m = %d, n1 = %d, n2 = %d, lt = %.2f, rd = %.2f" %\
                        (m, i, 2*m - i - 1, lt, rd))
                    print("left = %.2f, right = %.2f, u1 = %.4f, u2 = %.4f" %\
                        (left, right, u1, u2))
                    print("=" * 20)

                return u1, u2
    return 0, 0


def plot_transversal(relation, m_list, omega_list, precision):
    for m in m_list:
        for j in range(2 * m):
            # нарисуем границы прямоугольников
            left = j * pi / 2.0 / l1
            right = (j + 1) * pi / 2.0 / l1
            down = (2 * m - 1 - j) * pi / 2.0 / l2
            up = (2 * m - j) * pi / 2.0 / l2
            plt.plot([right, left, left, right, right],
                    [down, down, up, up, down], "k:")

            u1_list, u2_list = [], []
            margin = 1e-9
            u1 = left + margin
            while u1 < right:
                u2 = bisection(lambda x: relation(u1, x),
                        down + margin, up - margin, precision)
                if u2:
                    u1_list.append(u1)
                    u2_list.append(u2)
                u1 += precision
            plt.plot(u1_list, u2_list, "k-")

        for omega in omega_list:
            u1_list = np.linspace(0, m * pi / l1, 100)
            u2_list =\
                [((omega / c) ** 2 * (e2 * m2 - e1 * m1) + u1 ** 2) ** 0.5\
                for u1 in u1_list]
            plt.plot(u1_list, u2_list, "k--")
            u1, u2 = transversal_wavenumbers(relation, m, omega, 1e-7)
            plt.plot([u1], [u2], "ro")

    plt.xlabel(r"$u_1, cm^{-1}$")
    plt.ylabel(r"$u_2, cm^{-1}$")

    if DEBUG:
        plt.show()
    else:
        if relation is e_relation:
            name = "e"
        elif relation is m_relation:
            name = "m"
        else:
            name = "wtf"
        plt.savefig(name + ".png")
    plt.cla()


def longitudinal_wavenumber(relation, m, n, omega, precision):
    u1, u2 = transversal_wavenumbers(relation, m, omega, precision)
    if (u1 > 0):
        sqr_h = (omega / c) ** 2 * e1 * m1 - u1 ** 2 - (pi * n / b) ** 2
        if (sqr_h > 0):
            h = sqr_h ** 0.5
            return h
    return 0


def plot_longitudinal(relation, m_list, n_list, omega_list, precision):
    for m in m_list:
        for n in n_list:
            h_list = [longitudinal_wavenumber(relation, m, n, omega_list[0],
                precision)]
            o_list = [omega_list[0]]
            for omega in omega_list:
                h = longitudinal_wavenumber(relation, m, n, omega, precision)
                if (h == 0) and (h_list[-1] == 0):
                    h_list[-1] = h
                    o_list[-1] = omega
                elif h >= h_list[-1]:
                    h_list.append(h)
                    o_list.append(omega)
            plt.plot(o_list, h_list, "k-")
    plt.xlabel(r"$\omega, rad/s$")
    plt.ylabel(r"$h, cm^{-1}$")
    if DEBUG:
        plt.show()
    else:
        if relation is e_relation:
            name = "e_h"
        elif relation is m_relation:
            name = "m_h"
        else:
            name = "wtf"
        plt.savefig(name + ".png")
    plt.cla()


if __name__ == '__main__':
    plot_transversal(e_relation, [1,2,3,4], [2e10, 4e10, 6e10, 8e10], 1e-3)
    plot_transversal(m_relation, [1,2,3,4], [2e10, 4e10, 6e10, 8e10], 1e-3)
    plot_longitudinal(e_relation, [2,3,4], [1,2,3],
            np.linspace(2e10, 13e10, 200), 1e-4)
    plot_longitudinal(m_relation, [2,3,4], [0,1,2],
            np.linspace(2e10, 13e10, 200), 1e-4)

