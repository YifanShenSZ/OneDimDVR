import argparse

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("p", type=float, help="momentum expectation")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    ''' Initialize '''
    args = parse_args()
    pmean = args.p
    psigma = pmean / 20
    qsigma = 0.5 / psigma
    qmean = -10 - 5 * qsigma
    left = qmean - 5 * qsigma
    right = -left
    pmin = pmean - 5 * psigma
    pmax = pmean + 5 * psigma
    ''' Do the job '''
    with open("initial.in", 'w') as f:
        print("Coordinate expectation:", file=f)
        print(qmean, file=f)
        print("Momentum expectation:", file=f)
        print(pmean, file=f)
    with open("OneDimDVR-scatter.in", 'w') as f:
        print("Job type: (wavefunction, transmission)", file=f)
        print("transmission", file=f)
        print("Mass:", file=f)
        print(2000, file=f)
        print("Total propagation time:", file=f)
        print(1000000, file=f)
        print("Time step:", file=f)
        print(0.1, file=f)
        print("Output interval:", file=f)
        print(1, file=f)
        print("Left boundary:", file=f)
        print(left, file=f)
        print("Right boundary:", file=f)
        print(right, file=f)
        print("Grid spacing:", file=f)
        print(0.1, file=f)
        print("Minimum wave number to be absorbed:", file=f)
        print(pmin, file=f)
        print("Maximum wave number to be absorbed:", file=f)
        print(pmax, file=f)
