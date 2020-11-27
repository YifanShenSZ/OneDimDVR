'''
Animate density
'''

from typing import List
import argparse
import scipy.io
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as anm

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("state", type=str, nargs='+', help="The electronic state(s) of interest")
    parser.add_argument("-s","--save", action='store_true', help="Save animation")
    args = parser.parse_args()
    return args

# NGrids vector x
# NGrids by NLevels by NSnapshots 3rd-order tensor y
# y[i,j,k] is value on j-th level at x[i] position at k-th snapshot
# when there is only 1 level, y can alternatively be passed as a matrix
# default FPS = 25
def Animate2D(x, y,
title='', xlabel='q [a.u.]', ylabel='density',
speed=1.0, show=True, save=False, FileName='2D'):
    fig, axes = plt.subplots(y.shape[1], 1, squeeze=False)
    xmin = numpy.amin(x); xmax = numpy.amax(x); ymin = numpy.amin(y); ymax = numpy.amax(y)
    lines = []
    for j in range(y.shape[1]):
        axes[j][0].set_title(title); axes[j][0].set_xlabel(xlabel); axes[j][0].set_ylabel(ylabel)
        axes[j][0].set_xlim(xmin, xmax); axes[j][0].set_ylim(ymin, ymax)
        line, = axes[j][0].plot(x, y[:, j, 0])
        lines.append(line)
    def animate(i):
        for j in range(y.shape[1]): lines[j].set_ydata(y[:, j, i])
        return lines
    ani=anm.FuncAnimation(fig, animate, y.shape[2], interval = 40.0 / speed, blit = True)
    if(save): ani.save(FileName + ".gif")
    if(show): plt.show()

if __name__ == "__main__":
    args = parse_args()
    with scipy.io.FortranFile("grids.out", 'r') as f:
        grids = f.read_reals()
    NGrids = grids.shape[0]
    with scipy.io.FortranFile("density" + args.state[0] + ".out", 'r') as f:
        NSnapshots = 0
        try:
            while True:
                f.read_reals()
                NSnapshots += 1
        except:
            pass
    densities = numpy.empty((NGrids, len(args.state), NSnapshots))
    for i in range(len(args.state)):
        with scipy.io.FortranFile("density" + args.state[i] + ".out", 'r') as f:
            for j in range(NSnapshots):
                densities[:, i, j] = f.read_reals()
    Animate2D(grids, densities, save=args.save, FileName="density")