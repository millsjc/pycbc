import argparse
import h5py

parser = argparse.ArgumentParser()
parser.add_argument("--injfile", type=str, help="name of the input file")
parser.add_argument("--fixed-snr", type=float, help="fixed SNR value")
args = parser.parse_args()

with h5py.File(args.injfile, "r+") as f:
    fixed_snr = args.fixed_snr
    scale_factor = f["network_optimal_snr"][()] / fixed_snr
    for key in f.keys():
        if "distance" in key:
            f[key][()] = scale_factor * f[key][()]
        elif "snr" in key:
            f[key][()] =  f[key][()] / scale_factor
