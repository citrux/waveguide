from physics import lm_dispersion_relation, le_dispersion_relation


def plot_transversal(relation, m_list, omega_list, precision):
    rel = hw.e_relation if relation=="e" else hw.m_relation
    data = hw.transversal_data(rel, m_list, omega_list, precision)
    for x, y in data["curves"]:
        plt.plot(x, y, "k-")
    for x, y, u, v in data["frequencies"]:
        plt.plot(x, y, "k-", linewidth=0.5)
        plt.plot(u, v, "k--", linewidth=0.5)
    for x, y in data["borders"]:
        plt.plot(x, y, "wh", markersize=3)
    for x, y in data["solutions"]:
        plt.plot(x, y, "ko", markersize=3)
    rel = sw.e_relation if relation=="e" else sw.m_relation
    data = sw.transversal_data(rel, m_list, omega_list, precision)
    for x, y in data["curves"]:
        plt.plot([-t for t in x], y, "k-")
    for x, y, u, v in data["frequencies"]:
        plt.plot([-t for t in x], y, "k-", linewidth=0.5)
        plt.plot([-t for t in u], v, "k--", linewidth=0.5)
    for x, y in data["borders"]:
        plt.plot(-x, y, "wh", markersize=3)
    for x, y in data["solutions"]:
        plt.plot(-x, y, "ko", markersize=3)
    plt.xlabel(r"$u_1^2, cm^{-2}$")
    plt.ylabel(r"$u_2^2, cm^{-2}$")
    save_file(plt, relation + ".pdf")
    plt.cla()


def plot_longitudinal(m_list, n_list, omega_list, precision):
    for m, n in [(1,0), (1,1), (2,0), (1,2), (2,1), (3,0)]:
            h_list, o_list = [], []
            for omega in omega_list:
                h = hw.longitudinal_wavenumber(hw.m_relation, m, n, omega, precision)
                if h > 0:
                    h_list.append(h)
                    o_list.append(omega)
                if len(h_list) and (h == 0):
                    break;
            plt.plot(o_list, h_list, "b-")
            # ------------------------------
            h_list, o_list = [], []
            for omega in omega_list:
                h = sw.longitudinal_wavenumber(sw.m_relation, m, n, omega, precision)
                if h > 0:
                    h_list.append(h)
                    o_list.append(omega)
                if len(h_list) and (h == 0):
                    break;
            plt.plot(o_list, h_list, "b--")
            if h_list:
                plt.annotate(r"$\mu_{%d%d}$" % (m,n), (o_list[-1], h_list[-1]))
            # print("plot m %d %d" %(m,n))

            # ------------------------------
            # ------------------------------
    for m, n in [(0,1), (0,2), (1,1), (2,1)]:
            if (n == 0):
                continue
            print(n)
            h_list, o_list = [], []
            for omega in omega_list:
                h = hw.longitudinal_wavenumber(hw.e_relation, m, n, omega, precision)
                if h > 0:
                    h_list.append(h)
                    o_list.append(omega)
                if len(h_list) and (h == 0):
                    break;
            plt.plot(o_list, h_list, "k-")
            # ------------------------------
            h_list, o_list = [], []
            for omega in omega_list:
                h = sw.longitudinal_wavenumber(sw.e_relation, m, n, omega, precision)
                if h > 0:
                    h_list.append(h)
                    o_list.append(omega)
                if len(h_list) and (h == 0):
                    break;
            plt.plot(o_list, h_list, "k--")
            if h_list:
                plt.annotate(r"$\varepsilon_{%d%d}$" % (m,n), (o_list[-1], h_list[-1]))

            # plt.annotate
            # print("plot e %d %d" %(m,n))
    plt.xlabel(r"$\omega, rad/s$")
    plt.ylabel(r"$h, cm^{-1}$")
    save_file(plt, "dispersion.pdf")
    plt.cla()


def plot_transversal_field(relation, m, n, omega, border_orientation="|"):
    '''
    Строит картину поля в поперечном сечении
    border_orientation -- ориентация границы раздела, может принимать значения
    "|" и "-"
    '''
    # формируем набор точек для построения полей в областях
    Y1, X1 = np.mgrid[0:b:91j, c:a:106j]
    Y2, X2 = np.mgrid[0:b:91j, 0:c:46j]

    # определяем поперечные и продольные волновые числа
    sqr_u1, sqr_u2 = transversal_wavenumbers(relation, m, omega, 1e-3)
    h = ((omega/sol)**2 * e1 * m1 - sqr_u1 - (pi*n/b) ** 2) ** 0.5

    if sqr_u1 < 0:
        # гиперболические волны
        u1, u2 = np.sqrt(-sqr_u1), np.sqrt(sqr_u2)
        if relation is lm_dispersion_relation:
            family = "lm"
            # левая область (2): от 0 до с
            E2x = -(h ** 2 + (pi * n / b) ** 2) / u2 / h *\
                np.sinh (u1 * (c-a)) / np.sin(u2 * c) * \
                np.cos(u2*X2) * np.sin(pi * n / b * Y2)
            E2y = pi * n / h / b * np.sinh (u1 * (c-a)) / np.sin(u2 * c) *\
                np.sin(u2*X2) * np.cos(pi * n / b * Y2)
            E2xy = np.sqrt(E2x*E2x + E2y*E2y)

            H2x =  0 * X2
            H2y = -e0 * e2 * omega / u2 * np.sinh (u1 * (c-a)) / np.sin(u2 * c) * \
                np.cos(u2*X2) * np.sin(pi * n / b * Y2)
            H2xy = np.sqrt(H2x*H2x + H2y*H2y)

            # правая область (1): от c до a
            E1x = (h ** 2 + (pi * n / b) ** 2) / u1 / h *\
                np.cosh(u1*(X1 - a)) * np.sin(pi * n / b * Y1)
            E1y = pi * n / h / b * np.sinh(u1*(X1 - a)) * np.cos(pi * n / b * Y1)
            E1xy = np.sqrt(E1x*E1x + E1y*E1y)

            H1x =  0 * X1
            H1y = e0 * e1 * omega / u1 *\
                np.cosh(u1*(X1 - a)) * np.sin(pi * n / b * Y1)
            H1xy = np.sqrt(H1x*H1x + H1y*H1y)

        elif relation is le_dispersion_relation:
            family = "le"
            # левая область (2): от 0 до с
            E2x = 0 * X2
            E2y = -omega * m0 * m2 / u2 * np.cosh(u1 * (c-a)) / np.cos(u2 * c) *\
                np.sin(u2*X2) * np.cos(pi * n / b * Y2)
            E2xy = np.sqrt(E2x*E2x + E2y*E2y)

            H2x = (h ** 2 + (pi * n / b) ** 2) / u2 / h *\
                np.cosh(u1 * (c-a)) / np.cos(u2 * c) *\
                np.sin(u2*X2) * np.cos(pi * n / b * Y2)
            H2y = -pi * n / h / b * np.cosh(u1 * (c-a)) / np.cos(u2 * c) *\
                np.cos(u2*X2) * np.sin(pi * n / b * Y2)
            H2xy = np.sqrt(H2x*H2x + H2y*H2y)

            # правая область (1): от c до a
            E1x = 0 * X1
            E1y = -omega * m0 * m2 / u1 * np.sinh(u1*(X1 - a)) *\
                np.cos(pi * n / b * Y1)
            E1xy = np.sqrt(E1x*E1x + E1y*E1y)

            H1x =  (h ** 2 + (pi * n / b) ** 2) / u1 / h *\
                np.sinh(u1*(X1 - a)) * np.cos(pi * n / b * Y1)
            H1y = -pi * n / h / b * np.cosh(u1*(X1 - a)) * np.sin(pi * n / b * Y1)
            H1xy = np.sqrt(H1x*H1x + H1y*H1y)

    else:
        # гармонические волны
        u1, u2 = np.sqrt(sqr_u1), np.sqrt(sqr_u2)
        if relation is lm_dispersion_relation:
            family = "lm"
            # левая область (2): от 0 до с
            E2x = -(h ** 2 + (pi * n / b) ** 2) / u2 / h *\
                np.sin (u1 * (c-a)) / np.sin(u2 * c) * \
                np.cos(u2*X2) * np.sin(pi * n / b * Y2)
            E2y = pi * n / h / b * np.sin (u1 * (c-a)) / np.sin(u2 * c) *\
                np.sin(u2*X2) * np.cos(pi * n / b * Y2)
            E2xy = np.sqrt(E2x*E2x + E2y*E2y)

            H2x =  0 * X2
            H2y = -e0 * e2 * omega / u2 * np.sin (u1 * (c-a)) / np.sin(u2 * c) * \
                np.cos(u2*X2) * np.sin(pi * n / b * Y2)
            H2xy = np.sqrt(H2x*H2x + H2y*H2y)

            # правая область (1): от c до a
            E1x = -(h ** 2 + (pi * n / b) ** 2) / u1 / h *\
                np.cos(u1*(X1 - a)) * np.sin(pi * n / b * Y1)
            E1y = pi * n / h / b * np.sin(u1*(X1 - a)) * np.cos(pi * n / b * Y1)
            E1xy = np.sqrt(E1x*E1x + E1y*E1y)

            H1x =  0 * X1
            H1y = -e0 * e1 * omega / u1 *\
                np.cos(u1*(X1 - a)) * np.sin(pi * n / b * Y1)
            H1xy = np.sqrt(H1x*H1x + H1y*H1y)

        elif relation is le_dispersion_relation:
            family = "le"
            # левая область (2): от 0 до с
            E2x = 0 * X2
            E2y = -omega * m0 * m2 / u2 * np.cos(u1 * (c-a)) / np.cos(u2 * c) *\
                np.sin(u2*X2) * np.cos(pi * n / b * Y2)
            E2xy = np.sqrt(E2x*E2x + E2y*E2y)

            H2x =  (h ** 2 + (pi * n / b) ** 2) / u2 / h *\
                np.cos(u1 * (c-a)) / np.cos(u2 * c) *\
                np.sin(u2*X2) * np.cos(pi * n / b * Y2)
            H2y = -pi * n / h / b * np.cos(u1 * (c-a)) / np.cos(u2 * c) *\
                np.cos(u2*X2) * np.sin(pi * n / b * Y2)
            H2xy = np.sqrt(H2x*H2x + H2y*H2y)

            # правая область (1): от c до a
            E1x = 0 * X1
            E1y = -omega * m0 * m2 / u1 * np.sin(u1*(X1 - a)) * np.cos(pi * n / b * Y1)
            E1xy = np.sqrt(E1x*E1x + E1y*E1y)

            H1x =  (h ** 2 + (pi * n / b) ** 2) / u1 / h *\
                np.sin(u1*(X1 - a)) * np.cos(pi * n / b * Y1)
            H1y = -pi * n / h / b * np.cos(u1*(X1 - a)) * np.sin(pi * n / b * Y1)
            H1xy = np.sqrt(H1x*H1x + H1y*H1y)


    Emax = max(E1xy.max(), E2xy.max())
    Hmax = max(H1xy.max(), H2xy.max())

    if border_orientation == "|":
            lw = 1.5 * E2xy / Emax
            plt.streamplot(X2, Y2, E2x, E2y, density=.4, color="k", linewidth=lw)

            lw = 1.5 * H2xy / Hmax
            plt.streamplot(X2, Y2, H2x, H2y, density=.4, color="b", linewidth=lw)

            lw = 1.5 * E1xy / Emax
            plt.streamplot(X1, Y1, E1x, E1y, density=.7, color="k", linewidth=lw)

            lw = 1.5 * H1xy / Hmax
            plt.streamplot(X1, Y1, H1x, H1y, density=.7, color="b", linewidth=lw)

    elif border_orientation == "-":
            lw = 1.5 * E2xy / Emax
            plt.streamplot(Y2, X2, E2y, E2x, density=.4, color="k", linewidth=lw)

            lw = 1.5 * H2xy / Hmax
            plt.streamplot(Y2, X2, H2y, H2x, density=.4, color="b", linewidth=lw)

            lw = 1.5 * E1xy / Emax
            plt.streamplot(Y1, X1, E1y, E1x, density=.7, color="k", linewidth=lw)

            lw = 1.5 * H1xy / Hmax
            plt.streamplot(Y1, X1, H1y, H1x, density=.7, color="b", linewidth=lw)

    plt.tick_params(
        axis='both',          # changes apply to both-axises
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        left='off',
        right='off',
        labelbottom='off',
        labeltop='off',
        labelleft='off',
        labelright='off') # labels along the bottom edge are off
    plt.plot([c,c], [0,b], "k-")
    plt.title("$%s_{%d%d},\ %.1e\ rad/s$" %(family, m, n, omega))
    plt.savefig("field_%s_%d_%d_%.1e.pdf" %(family, m, n, omega))
    plt.cla()

