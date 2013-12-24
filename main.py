import os
import matplotlib.pyplot as plt
import numpy as np
import harmonic_waves as hw
import surface_waves as sw
from settings import *
from plots import plot_longitudinal


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=5)


def save_file(plot, fname):
    if os.path.isfile(fname):
        os.rename(fname, fname + ".old")
    plot.savefig(fname)


def wave_type(relation, m, omega):
    if relation is "lm":
        u2_max = np.pi * m / l2
    elif relation is "le":
        precision = 1e-4
        u2_max = hw.bisection(lambda x: np.tan(x*l2) + m1/m2*x*l1,
                np.pi * (2 * m - 1) / 2 / l2 + precision,
                np.pi * m / l2 - precision, precision)
    else:
        print("wtf?")
    if u2_max < omega / sol * (e2 * m2) ** .5:
        return sw
    else:
        return hw


def cutoff_omega(relation, m, n, omega_min=0, omega_max=1e11):
    def f(omega, relation=relation):
        wt = wave_type(relation, m, omega)
        if relation is "lm":
            relation = wt.e_relation
        elif relation is "le":
            relation = wt.m_relation
        u1, u2 = wt.transversal_wavenumbers(relation, m, omega, 1e-3)
        return (omega / sol) ** 2 * e2 * m2 - u2 ** 2 - (np.pi * n / b) ** 2
    return sw.bisection(f, omega_min, omega_max, 1e7)



if __name__ == '__main__':
    freqs = []
    for i in range(4):
        for j in range(4):
            freqs.append(("lm", i,j,cutoff_omega("lm", i, j)))
            freqs.append(("le", i,j,cutoff_omega("le", i, j)))

    for i in sorted(freqs, key=lambda x: x[3]):
        print("%s %d %d: %e" % i)
    # plot_transversal("e", [0,1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-3)
    # plot_transversal("m", [1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-3)
    plot = plot_longitudinal([0,1,2,3], [0,1,2], np.arange(0, 1e11, 1e8), 1e-3)
    save_file(plot, "dispersion.pdf")
