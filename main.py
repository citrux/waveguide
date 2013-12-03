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
        plt.plot(x, y, "k-")
    for x, y, u, v in data["frequencies"]:
        plt.plot(x, y, "k-", linewidth=0.5)
        plt.plot(u, v, "k--", linewidth=0.5)
    for x, y in data["borders"]:
        plt.plot(x, y, "wh")
    for x, y in data["solutions"]:
        plt.plot(x, y, "ko")
    rel = sw.e_relation if relation=="e" else sw.m_relation
    data = sw.transversal_data(rel, m_list, omega_list, precision)
    for x, y in data["curves"]:
        plt.plot([-t for t in x], y, "k-")
    for x, y, u, v in data["frequencies"]:
        plt.plot([-t for t in x], y, "k-", linewidth=0.5)
        plt.plot([-t for t in u], v, "k--", linewidth=0.5)
    for x, y in data["borders"]:
        plt.plot(-x, y, "wh")
    for x, y in data["solutions"]:
        plt.plot(-x, y, "ko")
    save_file(plt, relation + ".pdf")
    plt.cla()


if __name__ == '__main__':
    plot_transversal("e", [0,1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-3)
    plot_transversal("m", [1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-3)
