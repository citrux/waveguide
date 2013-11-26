# coding: utf-8
###############################################################################
# графическое решение дисперсионного уравнения
# совместимо и с python 2, и с python 3
# частоты СВЧ 3 ÷ 30 ГГц → omega = 20 ÷ 200 Град/с
###############################################################################

DEBUG = False

import timeit
import os
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

from math import pi, tan, log

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

a = 5.0
b = 3.0
c = 1.5
e1, m1, l1 = 1.0, 1.0, a - c
e2, m2, l2 = 5.0, 1.0, c
# скорость света, см/с
sol = 3e10
e0 = 8.85e-12
m0 = 4 * pi * 1e-7


def save_file(fname):
    if os.path.isfile(fname):
        os.rename(fname, fname + ".old")
    plt.savefig(fname)


def e_relation(u1, u2):
    return u1 * tan(u1 * l1) / e1 + u2 * tan(u2 * l2) / e2


def m_relation(u1, u2):
    return m1 * tan(u1 * l1) / u1 + m2 * tan(u2 * l2) / u2


def bisection(f, left, right, precision):
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
        return bisection(f, left, center, precision)
    return bisection(f, center, right, precision)


def transversal_wavenumbers(relation, m, omega, precision):
    for i in range(2*m):
        lt = ((2*m - i) / l2) ** 2 - (i / l1) ** 2 -\
                (2 * omega / sol / pi) ** 2 * (e2 * m2 - e1 * m1)
        rd = ((2*m - i - 1) / l2) ** 2 - ((i + 1) / l1) ** 2 -\
                (2 * omega / sol / pi) ** 2 * (e2 * m2 - e1 * m1)

        if lt * rd < 0:
            first = lambda x: x
            second = lambda x:\
                    (x ** 2 + (omega / sol) ** 2 * (e2 * m2 - e1 * m1)) ** 0.5

            left = i * pi / 2.0 / l1
            right = (i + 1) * pi / 2.0 / l1
            down = (2 * m - i - 1) * pi / 2.0 / l2
            up = (2 * m - i) * pi / 2.0 / l2

            new_left = down ** 2 -  (omega / sol) ** 2 * (e2 * m2 - e1 * m1)
            new_right = up ** 2 -  (omega / sol) ** 2 * (e2 * m2 - e1 * m1)
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
            # отмечаем решение
            u1, u2 = transversal_wavenumbers(relation, m, omega, precision)
            if u2:
                plt.plot([u1], [u2], "ko")

    for omega in omega_list:
        # рисуем гиперболу для частоты
        n = 0 if relation is m_relation else 1
        sqr_border = (omega/sol)**2 * e1 * m1 - (pi * n / b) ** 2
        border = sqr_border ** 0.5 if sqr_border > 0 else 0
        u1_list = np.arange(0, border, precision)
        u2_list =\
            [((omega / sol) ** 2 * (e2 * m2 - e1 * m1) + u1 ** 2) ** 0.5\
            for u1 in u1_list]
        plt.plot(u1_list, u2_list, "k-", linewidth=.5)
        u1_list = np.arange(0, m * pi / l1, precision)
        u2_list =\
            [((omega / sol) ** 2 * (e2 * m2 - e1 * m1) + u1 ** 2) ** 0.5\
            for u1 in u1_list]
        plt.plot(u1_list, u2_list, "k--", linewidth=.5)

        # рисуем отсечки для заданной частоты
        n_max = int(omega / sol * (e1 * m1) ** 0.5 * b / pi)
        while n <= n_max:
            sqr_border = (omega/sol)**2 * e1 * m1 - (pi * n / b) ** 2
            if sqr_border > 0:
                u1 = sqr_border ** 0.5
                u2 = ((omega / sol) ** 2 * (e2 * m2 - e1 * m1) +\
                        u1 ** 2) ** 0.5
                plt.plot([u1], [u2], "wh", linewidth=.5)
                plt.annotate('$n=%d$' % n, xy=(u1, u2+.05), ha="center",
                        va="bottom")
            n += 1

        # подписываем частоты
        src="%e" % omega
        mantissa, exponent = src.split("e")
        mantissa = mantissa.strip("0 .")
        exponent = exponent.strip(" +")
        plt.annotate(
                '$\\omega = %s \\cdot 10^{%s} rad/s$' % (mantissa, exponent),
                xy=(u1_list[-1]+.2, u2_list[-1]-.3),
                ha="right", va="top")

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
        save_file(name + ".pdf")
    plt.cla()


def longitudinal_wavenumber(relation, m, n, omega, precision):
    u1, u2 = transversal_wavenumbers(relation, m, omega, precision)
    if (u1 > 0):
        sqr_h = (omega / sol) ** 2 * e1 * m1 - u1 ** 2 - (pi * n / b) ** 2
        if (sqr_h > 0):
            h = sqr_h ** 0.5
            return h
    return 0


def plot_longitudinal(relation, m_list, n_list, omega_list, precision):
    if relation is e_relation:
        family = "varepsilon"
        name = "e_h"
        color = "black"
    elif relation is m_relation:
        family = "mu"
        name = "m_h"
        color = "blue"
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
            if h_list[-1] > 0:
                plt.plot(o_list, h_list, color=color)
                plt.annotate('$\\%s_{%d%d}$' % (family, m, n),
                        xy=(o_list[-1]+.1e10, h_list[-1]), ha="left",
                        va="center", color=color)
    plt.xlabel(r"$\omega, rad/s$")
    plt.ylabel(r"$h, cm^{-1}$")
    if DEBUG:
        plt.show()
    else:
        save_file(name + ".pdf")
    plt.cla()


def plot_transversal_field(relation, m, n, omega):
    u1, u2 = transversal_wavenumbers(relation, m, omega, 1e-3)
    h = ((omega/sol)**2 * e1 * m1 - u1 **2 - (pi*n/b) ** 2) ** 0.5
    sqr_g1 = u1 **2 + (pi * n / b) ** 2
    sqr_g2 = u2 **2 + (pi * n / b) ** 2
    if relation is e_relation:
        E1 = 1
        H1 = -E1 * e1 * e0 * omega * pi * n / u1 / h / b
        E2 = E1 * np.sin(u1 * (c - a)) / np.sin(u2 * c)
        H2 = H1 * np.cos(u1 * (c - a)) / np.cos(u2 * c)
        family = "e"
    elif relation is m_relation:
        H1 = 1
        E1 = H1 * m1 * m0 * omega * pi * n / u1 / h / b
        E2 = E1 * np.sin(u1 * (c - a)) / np.sin(u2 * c)
        H2 = H1 * np.cos(u1 * (c - a)) / np.cos(u2 * c)
        family = "m"

    # левая область (2): от 0 до с
    Y2, X2 = np.mgrid[0:b:91j, 0:c:46j]
    E2x = -(h * u2 * E2 - omega * m0 * m2 * pi * n / b * H2) * np.cos(u2*X2) *\
            np.sin(pi * n / b * Y2) / sqr_g2
    E2y = -(h * pi * n / b * E2 + omega * m0 * m2 * u2 * H2) * np.sin(u2*X2) *\
            np.cos(pi * n / b * Y2) / sqr_g2
    E2xy = np.sqrt(E2x*E2x + E2y*E2y)

    H2x = (omega * e0 * e2 * pi * n / b * E2 + h * u2 * H2) * np.sin(u2*X2) *\
            np.cos(pi * n / b * Y2) / sqr_g2
    H2y = -(omega * e0 * e2 * u2 * E2 - h * pi * n / b * H2) * np.cos(u2*X2) *\
            np.sin(pi * n / b * Y2) / sqr_g2
    H2xy = np.sqrt(H2x*H2x + H2y*H2y)

    # правая область (1): от c до a
    Y1, X1 = np.mgrid[0:b:91j, c:a:106j]
    E1x = -(h * u1 * E1 - omega * m0 * m1 * pi * n / b * H1) *\
            np.cos(u1*(X1 - a)) * np.sin(pi * n / b * Y1) / sqr_g1
    E1y = -(h * pi * n / b * E1 + omega * m0 * m1 * u1 * H1) *\
            np.sin(u1*(X1 - a)) * np.cos(pi * n / b * Y1) / sqr_g1
    E1xy = np.sqrt(E1x*E1x + E1y*E1y)

    H1x = (omega * e0 * e1 * pi * n / b * E1 + h * u1 * H1) *\
            np.sin(u1*(X1-a)) * np.cos(pi * n / b * Y1) / sqr_g1
    H1y = -(omega * e0 * e1 * u1 * E1 - h * pi * n / b * H1) *\
            np.cos(u1*(X1-a)) * np.sin(pi * n / b * Y1) / sqr_g1
    H1xy = np.sqrt(H1x*H1x + H1y*H1y)

    Emax = max(E1xy.max(), E2xy.max())
    Hmax = max(H1xy.max(), H2xy.max())

    lw = 1.5 * E2xy / Emax
    plt.streamplot(X2, Y2, E2x, E2y, density=.4, color='k', linewidth=lw)

    lw = 1.5 * H2xy / Hmax
    plt.streamplot(X2, Y2, H2x, H2y, density=.4, color='b', linewidth=lw)

    lw = 1.5 * E1xy / Emax
    plt.streamplot(X1, Y1, E1x, E1y, density=.7, color='k', linewidth=lw)

    lw = 1.5 * H1xy / Hmax
    plt.streamplot(X1, Y1, H1x, H1y, density=.7, color='b', linewidth=lw)

    plt.savefig("field_%s_%d_%d_%.1e.pdf" %(family, m, n, omega))
    plt.cla()


if __name__ == '__main__':

    start = timeit.default_timer()
    plot_transversal(e_relation, [1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-4)
    stop = timeit.default_timer()
    print(stop - start)
    start = timeit.default_timer()
    plot_transversal(m_relation, [1,2,3,4], [3e10, 5e10, 7e10, 9e10], 1e-4)
    stop = timeit.default_timer()
    print(stop - start)
    #plot_longitudinal(e_relation, [2,3], [1,2],
            #np.linspace(2e10, 10e10, 200), 1e-4)
    #plot_longitudinal(m_relation, [2,3], [0,1,2],
            #np.linspace(2e10, 10e10, 200), 1e-4)
    #plot_transversal_field(e_relation, 3, 1, 6e10)
