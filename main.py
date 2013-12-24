import os
import matplotlib.pyplot as plt
import numpy as np
import harmonic_waves as hw
import surface_waves as sw
from plots import plot_transversal, plot_longitudinal
from settings import *


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


if __name__ == '__main__':
    print(wave_type("lm", 2, 2e10))
    print(wave_type("lm", 2, 7e10))
    # plot_transversal("e", [0,1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-3)
    # plot_transversal("m", [1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-3)
    # plot_longitudinal([0,1,2,3], [0,1,2], np.arange(0, 1e11, 1e8), 1e-3)
