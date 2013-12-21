import os
import matplotlib.pyplot as plt
import numpy as np
import harmonic_waves as hw
import surface_wave as sw
from plots import plot_transversal, plot_longitudinal


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=5)


def save_file(plot, fname):
    if os.path.isfile(fname):
        os.rename(fname, fname + ".old")
    plot.savefig(fname)



if __name__ == '__main__':
    pass
    # plot_transversal("e", [0,1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-3)
    # plot_transversal("m", [1,2,3], [2e10, 4e10, 6e10, 8e10], 1e-3)
    # plot_longitudinal([0,1,2,3], [0,1,2], np.arange(0, 1e11, 1e8), 1e-3)
