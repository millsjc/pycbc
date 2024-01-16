'''Notes:
- get_indices_3_ifo_2 is the fastest and most memory efficient
BUT get_indices_3_ifo_3 is identical in speed and memory usage, but currently only does the middle parts. I thought this would be faster than it is.

TRY: https://pythonspeed.com/articles/speeding-up-numba/
'''
from contextlib import contextmanager
import numpy as np
from numba import njit
import time
import memory_profiler
from subprocess import Popen
from os import getpid
from time import sleep
from signal import SIGINT

from num_combinations import number_coincident_combinations_parts

# define a context manager to measure time and memory usage
@contextmanager
def measure_time_memory():
    t1_start = time.perf_counter()
    mem1_start = memory_profiler.memory_usage()[0]
    yield # this is where the body of the with statement goes, equivalent to what happens between the __enter__ and __exit__ methods of a class that is used with the with statement
    t1_stop = time.perf_counter()
    mem1_stop = memory_profiler.memory_usage()[0]
    print(
        "Elapsed time: %.2f [s]" % (t1_stop - t1_start)
    )
    print(
        "Memory used: %.2f [MiB]" % (mem1_stop - mem1_start)
    )

@contextmanager
def perf_stat():
    ''''Requires perf (Linux only) to be installed.'''
    try:
        p = Popen(["perf", "stat", "-p", str(getpid())])
        sleep(0.5)
        yield
        p.send_signal(SIGINT)
    except:
        pass
        

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

@njit
def product(a, b):
    return np.array([[i, j] for i in a for j in b])

@njit
def get_indices_3_ifo_3(tlen, t3_coinc_window, i_j_vals, n_tail, n_middle, dtype=np.int64):
    num_combinations = 2 * n_tail + n_middle
    idx = np.empty((num_combinations, 3), dtype=dtype)
    
    row_num = n_tail
    # £FIXME: Only does the middle part atm
    while row_num < num_combinations - n_tail:
        for i in range(t3_coinc_window + 1, tlen - t3_coinc_window - 1):
            # for j in prange(-t2_coinc_window, t2_coinc_window+1):
            #     for k in prange(-t2_coinc_window, t2_coinc_window+1):
            for j, k in i_j_vals:
                idx[row_num] = [i, j+i, k+i]
                row_num += 1
    return idx
# @njit
# def get_indices_3_ifo_3(tlen, t2_coinc_window, t3_coinc_window, n_tail, n_middle, dtype=np.int64):
#     num_combinations = 2 * n_tail + n_middle
#     idx = np.empty((num_combinations, 3), dtype=dtype)

#     i_vals = np.arange(t3_coinc_window + 1, tlen - t3_coinc_window - 1)[:, None, None]
#     j_vals = np.arange(-t2_coinc_window, t2_coinc_window + 1)[None, :, None]
#     k_vals = np.arange(-t3_coinc_window, t3_coinc_window + 1)[None, None, :]

#     # reshape i_vals for broadcasting
#     i_vals_reshaped = np.reshape(i_vals, (i_vals.shape[0], 1, 1))

#     indices = np.stack((i_vals_reshaped, j_vals + i_vals_reshaped, k_vals + i_vals_reshaped), axis=-1)
    
#     idx[n_tail:num_combinations-n_tail] = indices.reshape(-1, 3)

#     return idx


def get_parser():
    import argparse
    parser = argparse.ArgumentParser(description='Compute 3-IFO indices.')
    parser.add_argument('tlen', type=int, help='length of time series')
    parser.add_argument('t2_coinc_window', type=int, help='coincidence window for T2')
    parser.add_argument('t3_coinc_window', type=int, help='coincidence window for T3')
    # parser.add_argument('-n', '--num_combinations', type=int, help='number of index combinations')
    # parser.add_argument('--n_tail', type=int, help='number of index combinations')
    # parser.add_argument('--n_middle', type=int, help='number of index combinations')
    parser.add_argument('--dtype', type=str, default='int64', help='data type of the indices')
    return parser

def check_args(args):
    assert args.t2_coinc_window < args.t3_coinc_window, "t2_coinc_window must be less than t3_coinc_window"

def get_num_combinations(tlen, t2_coinc_window, t3_coinc_window):
    n_tail, n_middle = number_coincident_combinations_parts(tlen, t2_coinc_window, t3_coinc_window)
    num_combinations = 2 * n_tail + n_middle
    return n_tail,n_middle,num_combinations

if __name__ == "__main__":
    parser = get_args()
    args = parser.parse_args()
    check_args(args)
    
    dtype = np.dtype(args.dtype)
    
    n_tail, n_middle, num_combinations = get_num_combinations(args.tlen, args.t2_coinc_window, args.t3_coinc_window)
    print("Number of combinations: ", num_combinations)

    # Call all functions once to compile them
    t_compilation_start = time.perf_counter()
    i_j_vals = product(np.arange(-1, 1+1), np.arange(-1, 1+1))
    get_indices_3_ifo_1(4, 1, 1, dtype=dtype)
    get_indices_3_ifo_2(4, 1, 1, 26, dtype=dtype)
    get_indices_3_ifo_3(4, 1, i_j_vals, 13, 0, dtype=dtype)
    t_compilation_end = time.perf_counter()
    print("Time taken to compile the functions: ", t_compilation_end - t_compilation_start)

    print("Implementation 1, get_indices_3_ifo_1:")
    with measure_time_memory():
        indices_1 = get_indices_3_ifo_1(args.tlen, args.t2_coinc_window, args.t3_coinc_window, dtype=dtype)

    with perf_stat():
        indices_1 = get_indices_3_ifo_1(args.tlen, args.t2_coinc_window, args.t3_coinc_window, dtype=dtype)
        
    print("Implementation 2, get_indices_3_ifo_2:")
    with measure_time_memory():
        indices_2 = get_indices_3_ifo_2(args.tlen, args.t2_coinc_window, args.t3_coinc_window, num_combinations, dtype=dtype)
    
    with perf_stat():
        indices_2 = get_indices_3_ifo_2(args.tlen, args.t2_coinc_window, args.t3_coinc_window, num_combinations, dtype=dtype)
        
    i_j_vals = product(np.arange(-args.t2_coinc_window, args.t2_coinc_window+1), np.arange(-args.t3_coinc_window, args.t3_coinc_window+1))
    print("Implementation 3, get_indices_3_ifo_3:")
    with measure_time_memory():
        indices_3 = get_indices_3_ifo_3(args.tlen, args.t3_coinc_window, i_j_vals, n_tail, n_middle, dtype=dtype)
    
    with perf_stat():
        indices_3 = get_indices_3_ifo_3(args.tlen, args.t3_coinc_window, i_j_vals, n_tail, n_middle, dtype=dtype)
        
    # Check if the outputs of all implementations are the same

    if (indices_1 == indices_2).all() and (indices_1 == indices_3).all():
        print("The outputs of all implementations are the same, e.g. see the first 3 indices:")
    else:
        print("The outputs of the implementations are different")

    # Check the end of the arrays are something sensible
    print(indices_1[-3:])
    print(indices_2[-3:])
    print(indices_3[-3:])

