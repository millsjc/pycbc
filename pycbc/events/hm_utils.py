""" This module contains functions for calculating and manipulating higher harmonic
2-filter triggers.
"""
import logging
import numpy, numpy.random
import pycbc.waveform, pycbc.filter, pycbc.types, pycbc.psd, pycbc.fft, pycbc.conversions
import numpy as np
from numba import njit


# FIXME: remove the following three functions and replace with coincident search functions.
def get_coinc_indexes(idx_dict, time_delay_idx):
    """Return the indexes corresponding to coincident triggers

    Parameters
    ----------
    idx_dict: dict
        Dictionary of indexes of triggers above threshold in each
        detector
    time_delay_idx: dict
        Dictionary giving time delay index (time_delay*sample_rate) for
        each ifo

    Returns
    -------
    coinc_idx: list
        List of indexes for triggers in geocent time that appear in
        multiple detectors
    """
    coinc_list = np.array([], dtype=int)
    for ifo in idx_dict.keys():
        # Create list of indexes above threshold in single detector in geocent
        # time. Can then search for triggers that appear in multiple detectors
        # later.
        if len(idx_dict[ifo]) != 0:
            coinc_list = np.hstack(
                [coinc_list, idx_dict[ifo] - time_delay_idx[ifo]]
            )
    # Search through coinc_idx for repeated indexes. These must have been loud
    # in at least 2 detectors.
    counts = np.unique(coinc_list, return_counts=True)
    coinc_idx = counts[0][counts[1] > 1]
    return coinc_idx


def get_coinc_triggers(snrs, idx, t_delay_idx):
    """Returns the coincident triggers from the longer SNR timeseries

    Parameters
    ----------
    snrs: dict
        Dictionary of single detector SNR time series
    idx: list
        List of geocentric time indexes of coincident triggers
    t_delay_idx: dict
        Dictionary of indexes corresponding to light travel time from
        geocenter for each detector

    Returns
    -------
    coincs: dict
        Dictionary of coincident trigger SNRs in each detector
    """
    coincs = {ifo: snrs[ifo][idx + t_delay_idx[ifo]] for ifo in snrs}
    return coincs


def snr_2_filter(snr_dict_dom, snr_dict_sub, index, threshold, time_delay_idx):
    """Calculate the 2 filter SNR for all coincident triggers above
    threshold

    Parameters
    ----------
    snr_dict: dict
        Dictionary of individual detector SNRs
    index: list
        List of indexes (geocentric) for which to calculate coincident
        SNR
    threshold: float
        Coincident SNR threshold. Triggers below this are cut
    time_delay_idx: dict
        Dictionary of time delay from geocenter in indexes for each
        detector

    Returns
    -------
    rho_2_filt: numpy.ndarray
        Coincident 2 filter SNR for surviving triggers
    index: list
        The subset of input indexes corresponding to triggers that
        survive the cuts
    """
    # Restrict the snr timeseries to just the interesting points
    coinc_triggers_dom = get_coinc_triggers(snr_dict_dom, index, time_delay_idx)
    coinc_triggers_sub = get_coinc_triggers(snr_dict_sub, index, time_delay_idx)
    # Calculate the coincident snr

    snr_dom_array = np.array(
        [coinc_triggers_dom[ifo] for ifo in coinc_triggers_dom.keys()]
    )
    snr_sub_array = np.array(
        [coinc_triggers_sub[ifo] for ifo in coinc_triggers_sub.keys()]
    )
    rho_2_filt = abs(np.sqrt(np.sum(
        snr_dom_array * snr_dom_array.conj() + 
        snr_sub_array * snr_sub_array.conj(), axis=0)))
    # Apply threshold
    thresh_indexes = rho_2_filt > threshold
    index = index[thresh_indexes]
    # coinc_triggers = get_coinc_triggers(snr_dict, index, time_delay_idx)
    rho_2_filt = rho_2_filt[thresh_indexes]
    rho_2_filt_rss = rho_2_filt
    return rho_2_filt, rho_2_filt_rss, index

@njit 
def three_det_sum_jit(t1, t2, t3, t2_coinc_window, t3_coinc_window):
    N = len(t1)
    network_snr_sq = np.array([
        t1[i] + t2[j] + t3[k]
            for i in range(N)
                for j in range(max(i-t2_coinc_window, 0), min(N, i+t2_coinc_window+1))
                    for k in range(max(i-t3_coinc_window, 0), min(N, i+t3_coinc_window+1))
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

@njit
def get_indices_jit(N, t2_coinc_window, t3_coinc_window):
    idx = np.array([
        [i,j,k]
            for i in range(N)
                for j in range(max(i-t2_coinc_window, 0), min(N, i+t2_coinc_window+1))
                    for k in range(max(i-t3_coinc_window, 0), min(N, i+t3_coinc_window+1))
                ])
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

def number_coincident_combinations(N, t2_coinc_window, t3_coinc_window):
    """Assumes t3_coinc_window>t2_coinc_window, and N is the number of time samples."""
    n_tail = number_tail(t2_coinc_window, t3_coinc_window)
    n_middle = (N - 2*(t3_coinc_window+1)) * (2 * t2_coinc_window + 1) * (2 * t3_coinc_window + 1) 
    n_c = 2 * n_tail + n_middle
    return n_c

def get_i_j_k(coinc_idx, N, t2_coinc_window, t3_coinc_window, precalculated_idxs=None):
    """Get the indices i,j,k of the detector timeseries corresponding to 
    coincident index coinc_idx. Assumes t3_coinc_window>t2_coinc_window, 
    and N is the number of analyzed time samples."""
    if coinc_idx < 0:
        raise ValueError("indices cannot be negative.")
    n_c = number_coincident_combinations(N, t2_coinc_window, t3_coinc_window)
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
            i, j, k = precalculated_idxs[coinc_idx - (n_c+1)] + N - (2*t3_coinc_window + 2)
    return i, j, k

def return_coinc_indices(N, t2_coinc_window, t3_coinc_window):
    if 2*max(t2_coinc_window, t3_coinc_window) > N:
        raise NotImplementedError("analyzed time must be larger than "
            "the coincident window")
    # forgetting the tails at first
    largest_window = max(t2_coinc_window, t3_coinc_window)
    idx_1_mid = np.arange(N - 2*(t3_coinc_window+1)).repeat( \
        (2*t2_coinc_window+1)*(2*t3_coinc_window+1) \
    ) + largest_window + 1
    idx_2_mid = idx_1_mid + np.tile(
        np.arange(-t2_coinc_window, t2_coinc_window+1).repeat(2*t3_coinc_window+1), 
        (N - 2*(t3_coinc_window+1))
    )
    idx_3_mid = idx_1_mid + np.tile(
        np.arange(-t3_coinc_window, t3_coinc_window+1), 
        (N - 2*(t3_coinc_window+1)) * (2*t2_coinc_window+1)
    )
    idx_1_ends, idx_2_ends, idx_3_ends = get_indices_jit(
        2*(t3_coinc_window+1), t2_coinc_window, t3_coinc_window).T
    # now get the tails
    n_start = int(len(idx_1_ends) / 2)
    idx_1 = np.concatenate((idx_1_ends[:n_start], idx_1_mid, idx_1_ends[-n_start:] + N - (2*t3_coinc_window + 2)))
    idx_2 = np.concatenate((idx_2_ends[:n_start], idx_2_mid, idx_2_ends[-n_start:] + N - (2*t3_coinc_window + 2)))
    idx_3 = np.concatenate((idx_3_ends[:n_start], idx_3_mid, idx_3_ends[-n_start:] + N - (2*t3_coinc_window + 2)))
    return np.array([idx_1, idx_2, idx_3]).T

# below version is optimal when there aren't many triggers above the coinc threshold
# def three_det_sum_and_threshold(
#     t1, t2, t3, t2_coinc_window, t3_coinc_window, threshold, precalculated_idxs=None
# ):
#     N = len(t1)
#     if 2*max(t2_coinc_window, t3_coinc_window) > N:
#         raise NotImplementedError("analyzed time must be larger than "
#             "the coincident window")
#     # FIXME: have to convert to numpy so that numba can interpret. 
#     # Find another way if it's desirable to save memory.
#     t1 = t1.numpy()
#     t1 = np.real(t1 * t1.conj())
#     t2 = t2.numpy()
#     t2 = np.real(t1 * t1.conj())
#     t3 = t3.numpy()
#     t3 = np.real(t1 * t1.conj())

#     network_snr_sq = three_det_sum_jit(t1, t2, t3, t2_coinc_window, t3_coinc_window)
#     mask = network_snr_sq > threshold**2
#     n_above_thresh = sum(mask)
#     logging.info("{} ({}%) possible coincs are above threshold".format(
#         n_above_thresh, 100 * float(n_above_thresh) / len(mask)
#     ))
#     coinc_idx = np.nonzero(mask)[0]
#     # later can improve the below loop by vectorizing the get_i_j_k function.
#     det_idx = np.array([
#         get_i_j_k(_i, N, t2_coinc_window, t3_coinc_window, 
#             precalculated_idxs=precalculated_idxs) 
#         for _i in coinc_idx])
#     return np.sqrt(network_snr_sq[mask]), det_idx

def three_det_sum_and_threshold(
    t1, t2, t3, idx, threshold
):
    # FIXME: have to convert to numpy so that numba can interpret. 
    # Find another way if it's desirable to save memory.
    t1 = t1.numpy()
    t1 = np.real(t1 * t1.conj())
    t2 = t2.numpy()
    t2 = np.real(t1 * t1.conj())
    t3 = t3.numpy()
    t3 = np.real(t1 * t1.conj())

    network_snr_sq = three_det_sum_idx_jit(t1, t2, t3, idx)
    mask = network_snr_sq > threshold**2
    n_above_thresh = sum(mask)
    logging.info("{} ({}%) possible coincs are above threshold".format(
        n_above_thresh, 100 * float(n_above_thresh) / len(mask)
    ))
    return np.sqrt(network_snr_sq[mask]), idx[mask]