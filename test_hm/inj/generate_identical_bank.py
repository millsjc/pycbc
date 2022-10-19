import sys, h5py
import numpy as np

outfile=sys.argv[-1]

params = {}
with h5py.File(outfile, "r+") as f:
    for p in ["mass1", "mass2", "spin1z", "spin2z"]:
        params[p] = np.array(f[p])

with h5py.File("".join(outfile.split(".")[0]+"_bank.hdf"), "w") as f:
    for p,v in params.items():
        f[p] = v
    f["template_id"] = range(len(v))