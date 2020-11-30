'''
Visualize expectations

Output to table format for future advanced plots
'''

import scipy.io
import numpy
import matplotlib.pyplot as plt

class PolDef(object):
    def __init__(self, input_line):
        self.q = []
        self.p = []
        for term in input_line:
            coord = int(term)
            if coord > 0:
                self.q.append(coord)
            else:
                self.p.append(-coord)
    
    def title(self) -> str:
        outstr = ''
        for _ in self.q: outstr += 'q'
        for _ in self.p: outstr += 'p'
        return outstr

if __name__ == "__main__":
    with open("expectation.in", 'r') as f: lines = f.readlines()
    PolDefs = []
    for line in lines:
        strs = line.split()
        PolDefs.append(PolDef(strs[1:]))
    with scipy.io.FortranFile("snapshots.out", 'r') as f:
        snapshots = f.read_reals()
    NSnapshots = snapshots.shape[0]
    with scipy.io.FortranFile("expectation.out", 'r') as f:
        expectations = f.read_reals()
    NStates = int(expectations.shape[0] / len(PolDefs) / NSnapshots)
    expectations = expectations.reshape((len(PolDefs), NStates, NSnapshots), order='F')
    # Output to table format for future advanced plots
    for i in range(NStates):
        with open("expectation-"+str(i)+".txt", 'w') as f:
            print("time", end='\t', file=f)
            for pol in PolDefs: print(pol.title(), end='\t', file=f)
            print(file=f)
            for j in range(NSnapshots):
                print(snapshots[j], end='\t', file=f)
                for term in expectations[:,i,j]: print(term, end='\t', file=f)
                print(file=f)
    # Make a convenient plot
    for i in range(len(PolDefs)):
        fig, axes = plt.subplots(NStates, 1, squeeze=False)
        for j in range(NStates):
            axis = axes[NStates - 1 - j][0]
            axis.set_xlim(snapshots[0], snapshots[NSnapshots - 1])
            axis.set_ylim(numpy.amin(expectations[i, j, :]), numpy.amax(expectations[i, j, :]))
            axis.plot(snapshots, expectations[i, j, :])
        plt.show()