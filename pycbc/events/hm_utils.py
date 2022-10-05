""" This module contains functions for calculating and manipulating higher harmonic
2-filter triggers.
"""
import numpy as np
from numba import njit

def angle_between_detectors(tau_12, tau_13, tau_23):
    """Calculate the angle between ifo2 and ifo3 viewed from ifo1.
    
    Parameters
    ----------
    tau_12: float
        Light travel time between ifo1 and ifo2 (arbitrary units).
    tau_13: float
        Light travel time between ifo1 and ifo3 (arbitrary units).
    tau_23: float
        Light travel time between ifo2 and ifo3 (arbitrary units).

    Returns
    -------
    alpha_23: float
        Angle (in radians) between ifo2 and ifo3 viewed from ifo1.
    """
    b, a, c = tau_12, tau_13, tau_23
    alpha_23 = np.arccos((a**2 + b**2 - c**2) / (2*a*b))
    return alpha_23

def check_time_delay_is_physical(
    n2, n3, N2, N3, alpha_23, rounding=True):
    """Check if the delay time between two detector combinations
    is geometrically possible. The units for time are arbritrary 
    but must be consistent. However if using rounding units must be 
    time indices. 
    
    Parameters
    ----------
    n2: float or array
        Delay time between ifo1 and ifo2.
    n3: float or array
        Delay time between ifo1 and ifo2.
    N2: float
        Light travel time between ifo1 and ifo2.
    N3: float
        Light travel time between ifo1 and ifo3.
    alpha_23: float
        Angle (in radians) between ifo2 and ifo3 viewed from ifo1.
    rounding: bool
        If True will include the point just outside the geometrical 
        region in the physical region.
        
    Returns
    cond: bool, or boolean array
        False if delay times are unphysical. 
    """
    N3 = float(N3)
    p = N3 / N2 * np.cos(alpha_23)
    q = N3 / N2 * np.sin(alpha_23)
    if rounding:
        # below line ensures we take the dots just outside 
        # the ellipse to be safe and account for rounding
        # FIXME: revise below
        n2_safe = n2 - np.sign(n2)*0.1
        n3_safe = n3 - np.sign(n3)*0.1
    # Eq 3.4 for the ellipse from https://arxiv.org/pdf/gr-qc/9509042.pdf
    cond = n3_safe**2 + n2_safe**2*(N3/N2)**2 - 2*p*n2_safe*n3_safe < (q*N2)**2
    return cond

@njit 
def three_det_sum_jit(t1, t2, t3, t2_coinc_window, t3_coinc_window):
    tlen = len(t1)
    network_snr_sq = np.array([
        t1[i] + t2[j] + t3[k]
            for i in range(tlen)
                for j in range(max(i-t2_coinc_window, 0), min(tlen, i+t2_coinc_window+1))
                    for k in range(max(i-t3_coinc_window, 0), min(tlen, i+t3_coinc_window+1))
                ])
    return network_snr_sq

@njit
def three_det_sum_idx_jit(t1, t2, t3, idx):
    """Most efficient if you are using numba and have the indexes."""
    temp = np.array([
        t1[i] + t2[j] + t3[k]
            for i,j,k in idx
                ])
    return temp

def get_index_array_dtype(max_index):
    # reduce the size of the index array where possible
    bits_dtypes = [(8, np.uint8), (16, np.uint16), (32, np.uint32), (64, np.uint64)]
    for bits, dtype in bits_dtypes: 
        if 2**bits > max_index+1: break
    return dtype

@njit
def get_indices_jit(tlen, t2_coinc_window, t3_coinc_window, dtype=np.int64):
    idx = np.array([
        [i,j,k]
            for i in range(tlen)
                for j in range(max(i-t2_coinc_window, 0), min(tlen, i+t2_coinc_window+1))
                    for k in range(max(i-t3_coinc_window, 0), min(tlen, i+t3_coinc_window+1))
                ], dtype=dtype)
    return idx

def number_tail(t2_coinc_window, t3_coinc_window):
    n_tail = (
        sum(
            (t2_coinc_window + 1 + np.arange(0, t2_coinc_window+1)) * (t3_coinc_window + 1 + np.arange(0, t2_coinc_window+1))
        ) + \
        sum(
            (2* t2_coinc_window + 1) * (t3_coinc_window + 1 + np.arange(t2_coinc_window+1, t3_coinc_window+1))
        )
    )
    return n_tail

def number_coincident_combinations(tlen, t2_coinc_window, t3_coinc_window):
    """Assumes t3_coinc_window>t2_coinc_window, and tlen is the number of time samples."""
    n_tail = number_tail(t2_coinc_window, t3_coinc_window)
    n_middle = (tlen - 2*(t3_coinc_window+1)) * (2 * t2_coinc_window + 1) * (2 * t3_coinc_window + 1) 
    n_c = 2 * n_tail + n_middle
    return n_c

def get_i_j_k(coinc_idx, tlen, t2_coinc_window, t3_coinc_window, precalculated_idxs=None):
    """Get the indices i,j,k of the detector timeseries corresponding to 
    coincident index coinc_idx. Assumes t3_coinc_window>t2_coinc_window, 
    and tlen is the number of analyzed time samples."""
    if coinc_idx < 0:
        raise ValueError("indices cannot be negative.")
    n_c = number_coincident_combinations(tlen, t2_coinc_window, t3_coinc_window)
    n_tail = number_tail(t2_coinc_window, t3_coinc_window)
    largest_window = max(t2_coinc_window, t3_coinc_window)
    if (coinc_idx > n_tail-1)&(coinc_idx < n_c-n_tail):
        i, remainder = divmod(coinc_idx - n_tail, (2*t2_coinc_window+1)*(2*t3_coinc_window+1))
        i += largest_window + 1
        j, remainder_2 = divmod(remainder, (2*t3_coinc_window+1))
        j += i - (t2_coinc_window)
        k = remainder_2 + i - (t3_coinc_window)
    else:
        if precalculated_idxs is None:
            precalculated_idxs = get_indices_jit(2*(t3_coinc_window+1), t2_coinc_window, t3_coinc_window)
        
        if not (coinc_idx > (n_tail-1)):
            # at the beginning
            i, j, k = precalculated_idxs[coinc_idx]
        else:
            # at the end
            i, j, k = precalculated_idxs[coinc_idx - (n_c+1)] + tlen - (2*t3_coinc_window + 2)
    return i, j, k

def index_combinations(tlen, t2_coinc_window, t3_coinc_window, dtype=np.int64):
    """This messy function is equivalent to calling get_indices_jit, but is 
    generally faster."""
    if 2*max(t2_coinc_window, t3_coinc_window) > tlen:
        return get_indices_jit(tlen, t2_coinc_window, t3_coinc_window, dtype)
    # forgetting the tails at first
    largest_window = max(t2_coinc_window, t3_coinc_window)
    idx_1_mid = np.arange(tlen - 2*(t3_coinc_window+1), dtype=dtype).repeat( \
        (2*t2_coinc_window+1)*(2*t3_coinc_window+1) \
    ) + largest_window + 1
    idx_2_mid = idx_1_mid + np.tile(
        np.arange(-t2_coinc_window, t2_coinc_window+1, dtype=dtype).repeat(2*t3_coinc_window+1), 
        (tlen - 2*(t3_coinc_window+1))
    )
    idx_3_mid = idx_1_mid + np.tile(
        np.arange(-t3_coinc_window, t3_coinc_window+1, dtype=dtype), 
        (tlen - 2*(t3_coinc_window+1)) * (2*t2_coinc_window+1)
    )
    idx_1_ends, idx_2_ends, idx_3_ends = get_indices_jit(
        2*(t3_coinc_window+1), t2_coinc_window, t3_coinc_window, dtype=dtype).T
    # now get the tails
    n_start = int(len(idx_1_ends) / 2)
    idx_1 = np.concatenate((idx_1_ends[:n_start], idx_1_mid, idx_1_ends[-n_start:] + tlen - (2*t3_coinc_window + 2)))
    idx_2 = np.concatenate((idx_2_ends[:n_start], idx_2_mid, idx_2_ends[-n_start:] + tlen - (2*t3_coinc_window + 2)))
    idx_3 = np.concatenate((idx_3_ends[:n_start], idx_3_mid, idx_3_ends[-n_start:] + tlen - (2*t3_coinc_window + 2)))
    return np.array([idx_1, idx_2, idx_3]).T

def three_det_sum_and_threshold(t1, t2, t3, idx, threshold):
    dtype = get_index_array_dtype(np.max(idx))
    network_snr_sq = three_det_sum_idx_jit(abs(t1.numpy())**2, 
        abs(t2.numpy())**2, abs(t3.numpy())**2, idx.astype(dtype))
    mask = network_snr_sq > threshold**2
    return np.sqrt(network_snr_sq[mask]), idx[mask]

def inner_complex(a, b):
    """Assumes a 2D array of detector SNRs where zeroth index selects a detector."""
    return abs(np.sum(a * b.conjugate(), axis=0))

def snr_2_filter_and_threshold(snr_dom, snr_sub_perp, idx, threshold):
    nifos = len(snr_dom)
    dom = np.array([snr_dom[i][idx[:,i]] for i in range(nifos)])
    sub = np.array([snr_sub_perp[i][idx[:,i]] for i in range(nifos)])
    network_rho_dom = np.linalg.norm(dom, axis=0)
    network_rho_sub_perp = inner_complex(dom, sub) / network_rho_dom
    snr_2_filter = np.sqrt(network_rho_sub_perp**2 + network_rho_dom**2)
    mask = snr_2_filter > threshold
    return snr_2_filter[mask], idx[mask], mask

def maximal_coinc_in_ifo(snrs, det_idx, ifo_i=0):
    """Choose the maximum network snr for each time point in ifo_i 
    (defaults to the zeroth detector).
    Returns
    -------
    snrs
    det_idx
    i_max: numpy.array
        The indices that maximize the snr.
    """
    i_max = []
    j_start, j_end = 0, 0
    for c in np.unique(det_idx[:,ifo_i], return_counts=True)[1]:
        j_end += c
        i_max.append(np.argmax(snrs[j_start:j_end]) + j_start)
        j_start += c
    i_max = np.array(i_max)
    return snrs[i_max], det_idx[i_max], i_max