import os
import matplotlib.pyplot as plt
import numpy as np
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
    for m in m_list:
        for n in n_list:
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
            print("plot m %d %d" %(m,n))
            
            # ------------------------------
            # ------------------------------
            
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
            print("plot e %d %d" %(m,n))
    plt.xlabel(r"$\omega, rad/s$")
    plt.ylabel(r"$h, cm^{-1}$")
    save_file(plt, "dispersion.pdf")
    plt.cla()


if __name__ == '__main__':
    # plot_transversal("e", [0,1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-3)
    # plot_transversal("m", [1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-3)
    plot_longitudinal([0,1], [0,1], np.arange(0, 1e11, 1e8), 1e-3)
