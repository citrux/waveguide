import os
import matplotlib.pyplot as plt
import harmonic_waves as hw
import surface_wave as sw

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def save_file(plot, fname):
    if os.path.isfile(fname):
        os.rename(fname, fname + ".old")
    plot.savefig(fname)


def plot_transversal(relation, m_list, omega_list, precision):
    rel = hw.e_relation if relation=="e" else hw.m_relation
    data = hw.transversal_data(rel, m_list, omega_list, precision)
    for x, y in data["curves"]:
        plt.plot([t**2 for t in x], [t**2 for t in y], "k-")
    for x, y, u, v in data["frequencies"]:
        plt.plot([t**2 for t in x], [t**2 for t in y], "k-", linewidth=0.5)
        plt.plot([t**2 for t in u], [t**2 for t in v], "k--", linewidth=0.5)
    for x, y in data["borders"]:
        plt.plot(x**2, y**2, "wh", markersize=3)
    for x, y in data["solutions"]:
        plt.plot(x**2, y**2, "ko", markersize=3)
    rel = sw.e_relation if relation=="e" else sw.m_relation
    data = sw.transversal_data(rel, m_list, omega_list, precision)
    for x, y in data["curves"]:
        plt.plot([-t**2 for t in x], [t**2 for t in y], "k-")
    for x, y, u, v in data["frequencies"]:
        plt.plot([-t**2 for t in x], [t**2 for t in y], "k-", linewidth=0.5)
        plt.plot([-t**2 for t in u], [t**2 for t in v], "k--", linewidth=0.5)
    for x, y in data["borders"]:
        plt.plot(-x**2, y**2, "wh", markersize=3)
    for x, y in data["solutions"]:
        plt.plot(-x**2, y**2, "ko", markersize=3)
    plt.xlabel(r"$u_1^2, cm^{-2}$")
    plt.ylabel(r"$u_2^2, cm^{-2}$")
    save_file(plt, relation + ".pdf")
    plt.cla()


if __name__ == '__main__':
    plot_transversal("e", [0,1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-3)
    plot_transversal("m", [1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-3)
