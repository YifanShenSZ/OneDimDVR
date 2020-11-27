'''
Animate Wigner distribution
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
    parser.add_argument("state", type=str, help="The electronic state of interest")
    parser.add_argument("-s","--save", action='store_true', help="Save animation")
    args = parser.parse_args()
    return args

# take in list x and tolerance tol
# return the indice list of elements with absolute value > tol
def Pick_Significance(x:List, tol:float) -> List:
    indice = []
    for i in range(x.shape[0]):
        if abs(x[i]) > tol: indice.append(i)
    return indice

# Plot a 3D surface animation
# NMeshGrids order vectors x & y
# NMeshGrids by NSnapshots matrix z
# (x[i], y[i], z[i,j]) is the coordinate of i-th point at j-th snapshot
# default FPS = 25
def Animate3D(x, y, z,
title='', xlabel='q [a.u.]', ylabel='p [a.u.]', zlabel='density',
colormap='seismic',
speed=1.0, show=True, save=True, FileName='3D'):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(title); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_zlabel(zlabel)
    xmin = numpy.amin(x); xmax = numpy.amax(x); ax.set_xlim(xmin, xmax)
    ymin = numpy.amin(y); ymax = numpy.amax(y); ax.set_ylim(ymin, ymax)
    zmin = numpy.amin(z); zmax = numpy.amax(z); ax.set_zlim(zmin, zmax)
    tol = 0.01 * max(abs(zmax), abs(zmin))
    indice = Pick_Significance(z[:, 0], tol) # Neglect small values
    xplot = x[indice]; yplot = y[indice]; zplot = z[:, 0][indice]
    ax.plot_trisurf(xplot, yplot, zplot, cmap=colormap, vmin=zmin, vmax=zmax)
    def animate(i):
       ax.clear()
       ax.set_title(title); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_zlabel(zlabel)
       ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax); ax.set_zlim(zmin, zmax)
       indice = Pick_Significance(z[:, i], tol)
       xplot = x[indice]; yplot = y[indice]; zplot = z[:, i][indice]
       ax.plot_trisurf(xplot, yplot, zplot, cmap=colormap, vmin=zmin, vmax=zmax)
    ani = anm.FuncAnimation(fig, animate, frames=z.shape[1], interval=40.0/speed)
    if(save): ani.save(FileName + ".gif")
    if(show): plt.show()

if __name__ == "__main__":
    args = parse_args()
    with scipy.io.FortranFile("Wigner_q.out", 'r') as f:
        grids = f.read_reals()
    NGrids = grids.shape[0]
    with scipy.io.FortranFile("Wigner_p.out", 'r') as f:
        momenta = f.read_reals()
    NMomenta = momenta.shape[0]
    with scipy.io.FortranFile("Wigner" + args.state + ".out", 'r') as f:
        NSnapshots = 0
        try:
            while True:
                f.read_reals()
                NSnapshots += 1
        except:
            pass
    with scipy.io.FortranFile("Wigner" + args.state + ".out", 'r') as f:
        Wigner = numpy.empty((int(NGrids * NMomenta), NSnapshots))
        for i in range(NSnapshots):
            Wigner[:, i] = f.read_reals()
    q, p = numpy.meshgrid(grids, momenta)
    q = q.ravel()
    p = p.ravel()
    Animate3D(q, p, Wigner,
        title="Wigner distribution on state " + args.state,
        save=args.save, FileName="Wigner" + args.state)