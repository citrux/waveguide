import numpy as np
import matplotlib.pyplot as plt
from settings import *


def lm_dispersion_relation(sqr_u1, sqr_u2):
    '''
        Дисперсионное соотношение для LM-волн
        Если u_1^2 < 0, то волна экспонентцциально затухает в первой области
    '''
    if sqr_u1 < 0:
        u1, u2 = np.sqrt(-sqr_u1), np.sqrt(sqr_u2)
        return -u1 * np.tanh(u1 * l1) / e1 + u2 * np.tan(u2 * l2) / e2
    else:
        u1, u2 = np.sqrt(sqr_u1), np.sqrt(sqr_u2)
        return u1 * np.tan(u1 * l1) / e1 + u2 * np.tan(u2 * l2) / e2


def le_dispersion_relation(sqr_u1, sqr_u2):
    '''
        Дисперсионное соотношение для LE-волн
        Если u_1^2 < 0, то волна экспонентцциально затухает в первой области
    '''
    if sqr_u1 < 0:
        u1, u2 = np.sqrt(-sqr_u1), np.sqrt(sqr_u2)
        return m1 * np.tanh(u1 * l1) / u1 + m2 * np.tan(u2 * l2) / u2
    else:
        u1, u2 = np.sqrt(sqr_u1), np.sqrt(sqr_u2)
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
    '''
    u_1^2 |      |           m=2       |
          |      |                     |
          -------|---------------------|------
          |      |                     |
          |      |   u_1^2 = u_2^2-omega^2(e2m2 - e1m1)
          |      |         /           |
          |      |        /            |
          |      |       /             |  m=2
       1  |      |      /              |
          | m=1  |     /               |
          |      |    /                |
          |      |   /                 |
          |      |  /                  |
          |      | /                   |
          -------|/---------------------------
       0  |      /       m=1           |
          |     /|                     |
          -----/-|---------------------|------
          |  0/  |            1        |  u_2^2
          |  /   |                     |
          | /    |                     |
          |/     |                     |
          /      |                     |
          | m=0  |                     |
          |      |                     |
    '''
    delta = omega ** 2 * (e2*m2 - e1*m1) / sol ** 2 # sqr_u2 - sqr_u1
    n = 2 * m - 1  # сумма индексов квадратов с возможным решением для
                   # синусоидальных волн

    # границы квадратов
    border1 = lambda i: (np.pi * i / 2 / l1) ** 2
    border2 = lambda i: (np.pi * i / 2 / l2) ** 2

    if relation is lm_dispersion_relation:
        # sqr_top_harmonic_wavelength
        sthw = border2(n+1)

    elif relation is le_dispersion_relation:
        # top_harmonic_wavelength
        sthw = (bisection(lambda u2: m2 * np.tan(u2*l2) + m1*l1*u2,
                        np.sqrt(border2(n)) + precision,
                        np.sqrt(border2(n+1))-precision, precision))  ** 2

    # проверка на наличие синусоидальной волны
    if sthw > (omega / sol) ** 2 * (e2*m2 - e1*m1):
        # синусоидальная волна существует, так как точка пересечения с осью
        # лежит ниже максимально возможного значения

        # определимся с квадратом, в котором будем искать решение
        for i in range(n+1):
            if (border2(i) - border1(n+1-i) - omega**2 * (e2*m2-e1*m1)) *\
               (border2(i+1) - border1(n-i) - omega**2 * (e2*m2-e1*m1)) < 0:
                break
        # определились: sqr_u2 лежит в пределах от border2(i) до border2(i+1)
        # осталось бисекцией найти значение
        sqr_u2 = bisection(lambda x: relation(x - delta, x),
                        border2(i) + precision, border2(i+1) - precision,
                        precision)
        sqr_u1 = sqr_u2 - delta
    else:
        # волна затухает в поперечном сечении, так как точка пересечения с осью
        # лежит выше максимально возможного значения
        sqr_u2 = bisection(lambda x: relation(x - delta, x),
                        border2(n) + precision, sthw, precision)
        sqr_u1 = sqr_u2 - delta
    return sqr_u1, sqr_u2

print(transversal_wavenumbers(lm_dispersion_relation, 1, 4e10, 1e-4))
print(transversal_wavenumbers(lm_dispersion_relation, 1, 5e10, 1e-4))
print(transversal_wavenumbers(lm_dispersion_relation, 1, 6e10, 1e-4))


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
    sqr_u1, sqr_u2 = transversal_wavenumbers(relation, m, omega, precision)
    sqr_h = (omega / sol) ** 2 * e1 * m1 - sqr_u1 - (pi * n / b) ** 2
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

