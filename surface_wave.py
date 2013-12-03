# coding: utf-8
###############################################################################
# графическое решение дисперсионного уравнения
# совместимо и с python 2, и с python 3
# частоты СВЧ 3 ÷ 30 ГГц → omega = 20 ÷ 200 Град/с
###############################################################################

DEBUG = False

import numpy as np

from math import pi, tan, log, tanh
from settings import *

def e_relation(u1, u2):
    return -u1 * tanh(u1 * l1) / e1 + u2 * tan(u2 * l2) / e2


def m_relation(u1, u2):
    return m1 * tanh(u1 * l1) / u1 + m2 * tan(u2 * l2) / u2


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
    first = lambda x: x
    second = lambda x:\
            ((omega / sol) ** 2 * (e2 * m2 - e1 * m1) - x ** 2) ** 0.5
    j = 0 if relation is m_relation else 1
    down = (2 * m - 1 + j) * pi / 2.0 / l2
    up = (2 * m + j) * pi / 2.0 / l2

    sqr_right = (omega / sol) ** 2 * (e2 * m2 - e1 * m1) - down ** 2
    sqr_left = (omega / sol) ** 2 * (e2 * m2 - e1 * m1) - up ** 2

    left = sqr_left ** 0.5 if sqr_left > 0 else 0
    right = sqr_right ** 0.5 if sqr_right > 0 else 0

    margin = 1e-7
    if relation(left + margin, second(left + margin)) >= 0:
        u1 = bisection(lambda x: relation(first(x),\
            second(x)), left + margin, right - margin, precision)
        u2 = second(u1)

        if DEBUG:
            print("m = %d, down = %.2f, up = %.2f" %\
                (m, down, up))
            print("left = %.2f, right = %.2f, u1 = %.4f, u2 = %.4f" %\
                (left, right, u1, u2))
            print("=" * 20)

        return u1, u2
    return 0, 0


def transversal_data(relation, m_list, omega_list, precision):
    result = {
        "borders": [],
        "solutions": [],
        "curves": [],
        "frequencies": []
    }
    for m in m_list:
        j = 0 if relation is m_relation else 1
        down = (2 * m - 1 + j) * pi / 2.0 / l2
        up = (2 * m + j) * pi / 2.0 / l2
        left = 0
        right = (max(omega_list) / sol) * (e2 * m2 - e1 * m1) ** 0.5

        u1_list, u2_list = [], []
        margin = 1e-9
        u1 = left + margin
        while u1 < right:
            u2 = bisection(lambda x: relation(u1, x),
                    down + margin, up - margin, precision)
            if u2:
                u1_list.append(u1 ** 2)
                u2_list.append(u2 ** 2)
            u1 += precision
        result["curves"].append((u1_list, u2_list))

        for omega in omega_list:
            # отмечаем решение
            u1, u2 = transversal_wavenumbers(relation, m, omega, precision)
            if u1 and u2:
                result["solutions"].append([u1 ** 2, u2 ** 2])
                
    for omega in omega_list:
        # рисуем окружность для частоты
        n = 0 if relation is m_relation else 1
        border = (pi * n / b) ** 2 - (omega/sol)**2 * e1 * m1
        if border < 0:
            border = 0
        u = [0, border]
        v = [(omega / sol) ** 2 * (e2 * m2 - e1 * m1) - sqr_u1 for sqr_u1 in u]
        x = [border, (omega / sol) ** 2 * (e2 * m2 - e1 * m1)]
        y = [(omega / sol) ** 2 * (e2 * m2 - e1 * m1) - sqr_u1 for sqr_u1 in x]
        result["frequencies"].append((x, y, u, v))


        # рисуем отсечки для заданной частоты
        # n = 0
        # n_max = 3
        # while n <= n_max:
        #     sqr_border = (pi * n / b) ** 2 - (omega/sol) ** 2 * e1 * m1
        #     if sqr_border > 0 and\
        #         sqr_border < (omega / sol) ** 2 * (e2 * m2 - e1 * m1):
        #         u1 = sqr_border ** 0.5
        #         u2 = ((omega / sol) ** 2 * (e2 * m2 - e1 * m1) -\
        #                 u1 ** 2) ** 0.5
        #         result["borders"].append([u1, u2])
        #     n += 1
    return result



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

    # start = timeit.default_timer()
    # plot_transversal(e_relation, [1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-4)
    # stop = timeit.default_timer()
    # print(stop - start)
    # start = timeit.default_timer()
    plot_transversal(m_relation, [1,2,3], [3e10, 5e10, 7e10], 1e-3)
    # print(transversal_wavenumbers(m_relation, 1, 3e10, 1e-3))
    # stop = timeit.default_timer()
    # print(stop - start)
    #plot_longitudinal(e_relation, [2,3], [1,2],
            #np.linspace(2e10, 10e10, 200), 1e-4)
    #plot_longitudinal(m_relation, [2,3], [0,1,2],
            #np.linspace(2e10, 10e10, 200), 1e-4)
    #plot_transversal_field(e_relation, 3, 1, 6e10)
