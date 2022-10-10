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

def check_time_delay_is_physical(n2, n3, N2, N3, alpha_23):
    """Check if the delay time between two detector combinations
    is geometrically possible. The units for time are arbritrary 
    but must be consistent. 
    
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

    Returns
    -------
    cond: bool, or boolean array
        False if delay times are unphysical.

    Notes
    -----
    See p. 7 of https://arxiv.org/pdf/gr-qc/9509042.pdf
    """
    N3 = float(N3)
    p = N3 / N2 * np.cos(alpha_23)
    q = N3 / N2 * np.sin(alpha_23)
    # Eq 3.4 for the ellipse from https://arxiv.org/pdf/gr-qc/9509042.pdf
    cond = n3**2 + n2**2*(N3/N2)**2 - 2*p*n2*n3 < (q*N2)**2
    return cond

# below function not currently used anywhere.
def delay_bounds(n2, N2, N3, alpha_23):
    """Returns the bounds on possible delay times between detectors 1 
    and 3, given a delay between 1 and 2. The units for time are 
    arbritrary but must be consistent.

    Parameters
    ----------
    n2: float
        Delay time between ifo1 and ifo2.
    N2: float
        Light travel time between ifo1 and ifo2.
    N3: float
        Light travel time between ifo1 and ifo3.
    alpha_23: float
        Angle (in radians) between ifo2 and ifo3 viewed from ifo1.

    Returns
    -------
    min_delay: float
    max_delay: float

    Notes
    -----
    See p. 7 of https://arxiv.org/pdf/gr-qc/9509042.pdf
    """
    theta = np.arccos(n2 / N2)
    min_max = (N3 * np.cos(theta + alpha_23), N3 * np.cos(theta - alpha_23))
    return sorted(min_max)

@njit
def three_det_sum_idx_jit(t1, t2, t3, idx):
    """Most efficient if you are using numba and have the indexes."""
    temp = np.array([
        t1[i] + t2[j] + t3[k]
            for i,j,k in idx
                ])
    return temp

@njit
def two_det_sum_idx_jit(t1, t2, idx):
    return np.array([t1[i] + t2[j] for i,j in idx])

def get_index_array_dtype(max_index):
    # reduce the size of the index array where possible
    bits_dtypes = [(8, np.uint8), (16, np.uint16), (32, np.uint32), (64, np.uint64)]
    for bits, dtype in bits_dtypes: 
        if 2**bits > max_index+1: break
    return dtype

@njit
def get_indices_jit_3_ifo(tlen, t2_coinc_window, t3_coinc_window, dtype=np.int64):
    idx = np.array([
        [i,j,k]
            for i in range(tlen)
                for j in range(max(i-t2_coinc_window, 0), min(tlen, i+t2_coinc_window+1))
                    for k in range(max(i-t3_coinc_window, 0), min(tlen, i+t3_coinc_window+1))
                ], dtype=dtype)
    return idx

@njit
def get_indices_jit_2_ifo(tlen, t2_coinc_window, dtype=np.int64):
    idx = np.array([
        [i,j]
            for i in range(tlen)
                for j in range(max(i-t2_coinc_window, 0), min(tlen, i+t2_coinc_window+1))
                ], dtype=dtype)
    return idx

def index_combinations(tlen, t2_coinc_window, t3_coinc_window, dtype=np.int64):
    """For three detectors, this is equivalent to calling get_indices_jit_3_ifo, but is 
    generally faster. For two detectors this just calls get_indices_jit_2_ifo directly."""
    if t3_coinc_window is None:
        return get_indices_jit_2_ifo(tlen, t2_coinc_window, dtype)
    elif 2*max(t2_coinc_window, t3_coinc_window) > tlen:
        return get_indices_jit_3_ifo(tlen, t2_coinc_window, t3_coinc_window, dtype)
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
    idx_1_ends, idx_2_ends, idx_3_ends = get_indices_jit_3_ifo(
        2*(t3_coinc_window+1), t2_coinc_window, t3_coinc_window, dtype=dtype).T
    # now get the tails
    n_start = int(len(idx_1_ends) / 2)
    idx_1 = np.concatenate((idx_1_ends[:n_start], idx_1_mid, idx_1_ends[-n_start:] + tlen - (2*t3_coinc_window + 2)))
    idx_2 = np.concatenate((idx_2_ends[:n_start], idx_2_mid, idx_2_ends[-n_start:] + tlen - (2*t3_coinc_window + 2)))
    idx_3 = np.concatenate((idx_3_ends[:n_start], idx_3_mid, idx_3_ends[-n_start:] + tlen - (2*t3_coinc_window + 2)))
    return np.array([idx_1, idx_2, idx_3]).T

def detector_sum_and_threshold(snr_2_filt_rss, idx, threshold):
    """
    Parameters
    ----------
    snr_2_filt_rss: numpy.array2d
        zeroth index picks detector.
    idx: numpy.array2d
        The index combinations to sum over. Zeroth index picks detector.
    threshold: float

    Returns
    -------
    network_snr_2_filt_rss: numpy.array
        The network SNR for the events that survive the cut.
    idx: numpy.array2d
        The index combinations that survived the cut.
    """
    snr_2_filt_rss = abs(snr_2_filt_rss)**2
    nifos = len(snr_2_filt_rss)
    dtype = get_index_array_dtype(np.max(idx))
    if nifos == 3:
        network_snr_sq = three_det_sum_idx_jit(
            snr_2_filt_rss[0], snr_2_filt_rss[1], snr_2_filt_rss[2], idx.astype(dtype))
    elif nifos == 2:
        network_snr_sq = two_det_sum_idx_jit(
            snr_2_filt_rss[0], snr_2_filt_rss[1], idx.astype(dtype))
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

def maximal_coinc_in_ifo(snrs, time_idx):
    """Choose the maximum network snr for each time point.
    Returns
    -------
    i_max: numpy.array
        The indices that maximize the snr.
    """
    i_max = []
    j_start, j_end = 0, 0
    for c in np.unique(time_idx, return_counts=True)[1]:
        j_end += c
        i_max.append(np.argmax(snrs[j_start:j_end]) + j_start)
        j_start += c
    i_max = np.array(i_max)
    return i_max