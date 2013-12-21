import numpy as np
import matplotlib.pyplot as plt
from settings import *


def lm_dispersion_relation(u1_sqr, u2_sqr):
    '''
        Дисперсионное соотношение для LM-волн
        Если u_1^2 < 0, то волна экспонентцциально затухает в первой области
    '''
    if u1_sqr < 0:
        u1, u2 = sqrt(-u1_sqr), sqrt(u2_sqr)
        return -u1 * np.tanh(u1 * l1) / e1 + u2 * np.tan(u2 * l2) / e2
    else:
        u1, u2 = sqrt(u1_sqr), sqrt(u2_sqr)
        return u1 * np.tan(u1 * l1) / e1 + u2 * np.tan(u2 * l2) / e2


def le_dispersion_relation(u1_sqr, u2_sqr):
    '''
        Дисперсионное соотношение для LE-волн
        Если u_1^2 < 0, то волна экспонентцциально затухает в первой области
    '''
    if u1_sqr < 0:
        u1, u2 = sqrt(-u1_sqr), sqrt(u2_sqr)
        return m1 * np.tanh(u1 * l1) / u1 + m2 * np.tan(u2 * l2) / u2
    else:
        u1, u2 = sqrt(u1_sqr), sqrt(u2_sqr)
        return m1 * np.tan(u1 * l1) / u1 + m2 * np.tan(u2 * l2) / u2


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
    if (m == 0):
        return 0, 0
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


def transversal_data(relation, m_list, omega_list, precision):
    result = {
        "borders": [],
        "solutions": [],
        "curves": [],
        "frequencies": []
    }
    for m in m_list:
        for j in range(2 * m):
            # нарисуем границы прямоугольников
            left = j * pi / 2.0 / l1
            right = (j + 1) * pi / 2.0 / l1
            down = (2 * m - 1 - j) * pi / 2.0 / l2
            up = (2 * m - j) * pi / 2.0 / l2

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
            if u2:
                result["solutions"].append([u1 ** 2, u2 ** 2])

    for omega in omega_list:
        # рисуем гиперболу для частоты
        n = 0 if relation is le_dispersion_relation else 1
        border = (omega/sol)**2 * e1 * m1 - (pi * n / b) ** 2
        if border < 0:
            border = 0
        x = [0, border]
        y = [(omega / sol) ** 2 * (e2 * m2 - e1 * m1) + sqr_u1 for sqr_u1 in x]
        u = [border, (m * pi / l1) ** 2]
        v = [(omega / sol) ** 2 * (e2 * m2 - e1 * m1) + sqr_u1 for sqr_u1 in u]
        result["frequencies"].append((x, y, u, v))

        # рисуем отсечки для заданной частоты
        n_max = int(omega / sol * (e1 * m1) ** 0.5 * b / pi)
        while n <= n_max:
            border = (omega/sol)**2 * e1 * m1 - (pi * n / b) ** 2
            if border > 0:
                sqr_u1 = border
                sqr_u2 = (omega / sol) ** 2 * (e2 * m2 - e1 * m1) + sqr_u1
                result["borders"].append([sqr_u1, sqr_u2])
            n += 1
    return result


def longitudinal_wavenumber(relation, m, n, omega, precision):
    u1, u2 = transversal_wavenumbers(relation, m, omega, precision)
    if (u1 > 0):
        sqr_h = (omega / sol) ** 2 * e1 * m1 - u1 ** 2 - (pi * n / b) ** 2
        if (sqr_h > 0):
            h = sqr_h ** 0.5
            return h
    return 0


def longitudinal_data(relation, m_list, n_list, omega_list, precision):
    result = {
        "curves": []
    }
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
            result["curves"].append((o_list, h_list));
    return result

