import sys, h5py

outfile=sys.argv[-1]

with h5py.File(outfile, "r+") as f:
    for p in ["spin1x", "spin2x", "spin1y", "spin2y"]:
        for i in range(len(f[p])):
            f[p][i] = 0