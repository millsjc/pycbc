import ctypes
import numpy as np
from numba import njit, prange
import argparse
import time
import memory_profiler

parser = argparse.ArgumentParser(description='Compute 3-IFO indices.')
parser.add_argument('tlen', type=int, help='length of time series')
parser.add_argument('t2_coinc_window', type=int, help='coincidence window for T2')
parser.add_argument('t3_coinc_window', type=int, help='coincidence window for T3')
parser.add_argument('num_combinations', type=int, help='number of index combinations')
args = parser.parse_args()

@njit
def get_indices_3_ifo_1(tlen, t2_coinc_window, t3_coinc_window, dtype=np.int64):
    idx = np.array(
        [
            [i, j, k]
            for i in range(tlen)
            for j in range(
                max(i - t2_coinc_window, 0), min(tlen, i + t2_coinc_window + 1)
            )
            for k in range(
                max(i - t3_coinc_window, 0), min(tlen, i + t3_coinc_window + 1)
            )
        ],
        dtype=dtype,
    )
    return idx

@njit
def get_indices_3_ifo_2(tlen, t2_coinc_window, t3_coinc_window, num_combinations, dtype=np.int64):
    idx = np.empty((num_combinations, 3), dtype=dtype)
    row_num = 0
    for i in range(tlen):
        for j in range(
            max(i - t2_coinc_window, 0), min(tlen, i + t2_coinc_window + 1)
        ):
            for k in range(
                max(i - t3_coinc_window, 0), min(tlen, i + t3_coinc_window + 1)
            ):
                idx[row_num] = [i, j, k]
                row_num += 1
    return idx


def get_indices_3_ifo_3(tlen, t2_coinc_window, t3_coinc_window, num_combinations, dtype=np.int64):
    # # Load the shared library containing the C function
    # lib = ctypes.cdll.LoadLibrary('./indices.so')
    # idx = np.empty((num_combinations, 3), dtype=dtype)

    # # Call the C function
    # get_indices = lib.get_indices_3_ifo_forloop
    # get_indices.argtypes = [ctypes.c_int64, ctypes.c_int64, ctypes.c_int64, np.ctypeslib.ndpointer(dtype=np.int64, shape=(num_combinations, 3)), ctypes.c_int64]
    # get_indices.restype = None
    # get_indices(tlen, t2_coinc_window, t3_coinc_window, idx, num_combinations)
    
    idx = get_indices_3_ifo_1(tlen, t2_coinc_window, t3_coinc_window, dtype=np.int64)
    
    return idx

@njit(parallel=True)
def get_indices_3_ifo_2(tlen, t2_coinc_window, t3_coinc_window, num_combinations, dtype=np.int64):
    idx = np.empty((num_combinations, 3), dtype=dtype)
    row_num = 0
    for i in prange(tlen):
        for j in prange(
            max(i - t2_coinc_window, 0), min(tlen, i + t2_coinc_window + 1)
        ):
            for k in prange(
                max(i - t3_coinc_window, 0), min(tlen, i + t3_coinc_window + 1)
            ):
                idx[row_num] = [i, j, k]
                row_num += 1
    return idx

# Measure the time and memory usage of the first implementation
t1_start = time.perf_counter()
mem1_start = memory_profiler.memory_usage()[0]
indices_1 = get_indices_3_ifo_1(args.tlen, args.t2_coinc_window, args.t3_coinc_window)
mem1_end = memory_profiler.memory_usage()[0]
t1_end = time.perf_counter()

# Measure the time and memory usage of the second implementation
t2_start = time.perf_counter()
mem2_start = memory_profiler.memory_usage()[0]
indices_2 = get_indices_3_ifo_2(args.tlen, args.t2_coinc_window, args.t3_coinc_window, args.num_combinations)
mem2_end = memory_profiler.memory_usage()[0]
t2_end = time.perf_counter()

# Measure the time and memory usage of the third implementation
t3_start = time.perf_counter()
mem3_start = memory_profiler.memory_usage()[0]
indices_3 = get_indices_3_ifo_3(args.tlen, args.t2_coinc_window, args.t3_coinc_window, args.num_combinations)
mem3_end = memory_profiler.memory_usage()[0]
t3_end = time.perf_counter()

# Check if the outputs of all implementations are the same

if (indices_1 == indices_2)&(indices_2 == indices_3).all():
    print("The outputs of all implementations are the same, e.g. see the first 3 indices:")
else:
    print("The outputs of the implementations are different")

#check the end of the arrays are something sensible
print(indices_1[-3:])
print(indices_2[-3:])
print(indices_3[-3:])

# Print the time and memory usage of both implementations
print("Time taken by implementation 1: ", t1_end - t1_start)
print("Memory used by implementation 1: ", mem1_end - mem1_start)
print("Time taken by implementation 2: ", t2_end - t2_start)
print("Memory used by implementation 2: ", mem2_end - mem2_start)
print("Time taken by implementation 3: ", t3_end - t3_start)
print("Memory used by implementation 3: ", mem3_end - mem3_start)

