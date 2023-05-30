import numpy as np
from numba import njit

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


# call the function with the following parameters
# get_indices_3_ifo_2(100, 8, 16, 50084)
# create a function that will display the (3, num_combinations) array 
# as three scrollable rows.
# Then iterate slowly through the rows visually highlighting the table sqaure in colour
# in the order sugegsted by the above function. The colours 
# should be unique for every unique difference in number of indices.