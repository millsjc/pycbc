"""Here we experiment with loading from memory instead of computing the indices on the fly
As the number of snr series that have to be anluysed increases, it matters less which indices method you use. 
For one of two templates, it is probably best to use the numpy / memmap method if speed is the issue, or joblib/h5py if memory is the issue.

Tested with:
    ntemplates = 10
    nchunks = 3
Results:
    numpy:
    Elapsed time: 0.15 [s]
    Memory used: 0.00 [MiB]
    Summing and thresholding
    Elapsed time: 67.24 [s]
    Memory used: 2064.80 [MiB]


    memmap:
    Elapsed time: 0.12 [s]
    Memory used: 0.01 [MiB]
    Summing and thresholding
    Elapsed time: 66.17 [s]
    Memory used: 490.28 [MiB]


    hdf5:
    Elapsed time: 2.47 [s]
    Memory used: -462.25 [MiB]
    Summing and thresholding
    Elapsed time: 66.32 [s]
    Memory used: 0.00 [MiB]


    joblib:
    Elapsed time: 2.16 [s]
    Memory used: 0.05 [MiB]
    Summing and thresholding
    Elapsed time: 66.64 [s]
    Memory used: 429.21 [MiB]
"""
import time
import numpy as np
from benchmark2 import get_parser, check_args, get_indices_3_ifo_2, get_num_combinations, measure_time_memory, perf_stat
import h5py
from joblib import dump, load
from pycbc.events import hm_utils

sum_and_threshold = hm_utils.detector_sum_and_threshold



def save_indices(indices, dt):
    # save as various types
    np.save('indices.npy', indices)
    
    mmap = np.memmap('indices_mmap.dat', dtype=indices.dtype, mode='w+', shape=indices.shape)
    
    mmap[:] = indices[:]


    with h5py.File('indices.hdf5', 'w') as f:
        f.create_dataset("indices", data=indices, compression="gzip", dtype=dt)

    dump(indices, 'indices.joblib', compress=True)
    
    
def load_memmap():
    mmap = np.memmap('indices_mmap.dat', dtype=indices.dtype, mode='r', shape=indices.shape)
    return mmap  # Return the memory-mapped array directly

def load_numpy():
    return np.load('indices.npy')

def load_hdf5():
    with h5py.File('indices.hdf5', 'r') as f:
        return f['indices'][:]
    
def load_joblib():
    return load('indices.joblib')

if __name__ == "__main__":
    args = get_parser().parse_args()
    check_args(args)
    
    dt = args.dtype if args.dtype is not None else "uint16"
    dtype = np.dtype(dt)
    
    snrs = 2.5+np.random.chisquare(2, size=(3, args.tlen))

    print("dtype: ", dtype)
    n_tail, n_middle, num_combinations = get_num_combinations(args.tlen, args.t2_coinc_window, args.t3_coinc_window)
    print("Number of combinations: ", num_combinations)

    # Call all functions once to compile them
    t_compilation_start = time.perf_counter()
    get_indices_3_ifo_2(4, 1, 1, 26, dtype=dtype)
    t_compilation_end = time.perf_counter()
    print("Time taken to compile the functions: ", t_compilation_end - t_compilation_start)

    print("Implementation 2, get_indices_3_ifo_2:")
    with measure_time_memory():
        indices = get_indices_3_ifo_2(args.tlen, args.t2_coinc_window, args.t3_coinc_window, num_combinations, dtype=dtype)
    
    print("Summing and thresholding")
    with measure_time_memory():
        net_snr, det_idx = sum_and_threshold(snrs, indices, 7)
    
    print(indices)
    
    save_indices(indices, dt)


    print("Loading from memory and doing addition")
    ntemplates = 10
    nchunks = 3
    print("\n")
    print("numpy:")
    with measure_time_memory():
        indices = load_numpy()
    print("Summing and thresholding")
    with measure_time_memory():
        for i in range(ntemplates * nchunks):
            snrs = 2.5+np.random.chisquare(2, size=(3, args.tlen))
            net_snr, det_idx = sum_and_threshold(snrs, indices, 7)

    print("\n")
    print("memmap:")
    with measure_time_memory():
        indices = load_memmap()
    print("Summing and thresholding")
    with measure_time_memory():
        for i in range(ntemplates * nchunks):
            snrs = 2.5+np.random.chisquare(2, size=(3, args.tlen))
            net_snr, det_idx = sum_and_threshold(snrs, indices, 7)

    print("\n")
    print("hdf5:")  
    with measure_time_memory():
        indices = load_hdf5()
    print("Summing and thresholding")
    with measure_time_memory():
        for i in range(ntemplates * nchunks):
            snrs = 2.5+np.random.chisquare(2, size=(3, args.tlen))
            net_snr, det_idx = sum_and_threshold(snrs, indices, 7)

    print("\n")
    print("joblib:")
    with measure_time_memory():
        indices = load_joblib()
    print("Summing and thresholding")
    with measure_time_memory():
        for i in range(ntemplates * nchunks):
            snrs = 2.5+np.random.chisquare(2, size=(3, args.tlen))
            net_snr, det_idx = sum_and_threshold(snrs, indices, 7)