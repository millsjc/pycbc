""" This module contains functions for calculating and manipulating higher
harmonic 2-filter triggers.
"""
from math import ceil

import numpy as np
from numba import njit

from pycbc import detector


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
    alpha_23 = np.arccos((a**2 + b**2 - c**2) / (2 * a * b))
    return alpha_23


def get_coincident_window(ifo_list, timing_error, sample_rate):
    nifo = len(ifo_list)
    if nifo > 3:
        raise NotImplementedError("Cannot have more than 3 detectors.")
    t2_coinc_window = None
    t3_coinc_window = None
    t23_coinc_window = None
    alpha_23 = None
    segment_overlap = 0
    if nifo > 1:
        # calculate coincident window.
        ifo1 = detector.Detector(ifo_list[0])
        ifo2 = detector.Detector(ifo_list[1])
        time_diff_12 = ifo1.light_travel_time_to_detector(ifo2)
        t2_coinc_window = int(ceil((time_diff_12 + timing_error) 
                                   * sample_rate))
        segment_overlap = t2_coinc_window
        if nifo > 2:
            ifo3 = detector.Detector(ifo_list[2])
            time_diff_13 = ifo1.light_travel_time_to_detector(ifo3)
            time_diff_23 = ifo2.light_travel_time_to_detector(ifo3)
            if time_diff_12 > time_diff_13:
                # FIXME: consider switching labels, ifo_list order and 
                # anything else necessary to force this condition
                raise NotImplementedError("t3_coinc_window must be "
                                          "larger than t2_coinc window.")
            alpha_23 = angle_between_detectors(time_diff_12, time_diff_13,
                                               time_diff_23)
            t3_coinc_window = int(ceil((time_diff_13 + timing_error) 
                                       * sample_rate))
            t23_coinc_window = int(ceil((time_diff_23 + timing_error) 
                                        * sample_rate))
            segment_overlap = max(segment_overlap, t3_coinc_window)
    return (t2_coinc_window, t3_coinc_window, segment_overlap, 
            alpha_23, t23_coinc_window)


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
    cond = n3**2 + n2**2 * (N3 / N2) ** 2 - 2 * p * n2 * n3 < (q * N2) ** 2
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


# following three functions aren't really used.
def number_tail(t2_coinc_window, t3_coinc_window):
    n_tail = sum(
        (t2_coinc_window + 1 + np.arange(0, t2_coinc_window + 1))
        * (t3_coinc_window + 1 + np.arange(0, t2_coinc_window + 1))
    ) + sum(
        (2 * t2_coinc_window + 1)
        * (t3_coinc_window + 1 + np.arange(t2_coinc_window + 1, 
                                           t3_coinc_window + 1))
    )
    return n_tail


def number_coincident_combinations(tlen, t2_coinc_window, t3_coinc_window):
    """Assumes t3_coinc_window>t2_coinc_window, and tlen is the number 
    of time samples.
    """
    n_tail = number_tail(t2_coinc_window, t3_coinc_window)
    T_2 = 2 * t2_coinc_window + 1
    T_3 = 2 * t3_coinc_window + 1
    n_middle = (tlen - 2 * (t3_coinc_window + 1)) * T_2 * T_3
    n_c = 2 * n_tail + n_middle
    return n_c


def get_i_j_k(
    coinc_idx, tlen, t2_coinc_window, t3_coinc_window, precalculated_idxs=None
):
    """Get the indices i,j,k of the detector timeseries corresponding to
    coincident index coinc_idx. Assumes 3 detectors, t3_coinc_window 
    > t2_coinc_window, and tlen is the number of analyzed time samples.
    """
    if coinc_idx < 0:
        raise ValueError("indices cannot be negative.")
    n_c = number_coincident_combinations(tlen, t2_coinc_window, 
                                         t3_coinc_window)
    n_tail = number_tail(t2_coinc_window, t3_coinc_window)
    largest_window = max(t2_coinc_window, t3_coinc_window)
    T_2 = 2 * t2_coinc_window + 1
    T_3 = 2 * t3_coinc_window + 1
    if (coinc_idx > n_tail - 1) & (coinc_idx < n_c - n_tail):
        i, remainder = divmod(coinc_idx - n_tail, T_2 * T_3)
        i += largest_window + 1
        j, remainder_2 = divmod(remainder, T_3)
        j += i - (t2_coinc_window)
        k = remainder_2 + i - (t3_coinc_window)
    else:
        if precalculated_idxs is None:
            precalculated_idxs = get_indices_jit_3_ifo(
                2 * (t3_coinc_window + 1), t2_coinc_window, t3_coinc_window
            )

        if not coinc_idx > (n_tail - 1):
            # at the beginning
            i, j, k = precalculated_idxs[coinc_idx]
        else:
            # at the end
            i, j, k = (
                precalculated_idxs[coinc_idx - (n_c + 1)]
                + tlen
                - (2 * t3_coinc_window + 2)
            )
    return i, j, k


def max_number_timeslides(
    idx_onsource, t2_coinc_window, t3_coinc_window, t23_coinc_window
):
    nifo = np.shape(idx_onsource)[-1]
    T_2 = 2 * t2_coinc_window + 1
    if nifo == 2:
        n_incr = T_2**2
    elif nifo == 3:
        T_3 = 2 * t3_coinc_window + 1
        T_23 = 2 * t23_coinc_window + 1
        n_incr = T_2 * T_3 * (max(T_2, T_3) + T_23)
    n_c = len(idx_onsource)
    n_slides = n_c / float(n_incr) - 1
    return int(n_slides)


def perform_timeslide(idx_onsource, slide_no, t2_coinc_window, 
                      t3_coinc_window, t23_coinc_window):
    """Returns a copy of idx_combinations with the time indexes slid 
    so they are no longer coincident"""
    idx = idx_onsource.copy()
    nifo = np.shape(idx)[-1]
    T_2 = 2 * t2_coinc_window + 1
    if nifo == 2:
        n2_incr = T_2**2 * slide_no
        assert n2_incr < len(idx), "too many timeslides"
        idx[:, 1] = np.roll(idx[:, 1], n2_incr)
    elif nifo == 3:
        T_3 = 2 * t3_coinc_window + 1
        T_23 = 2 * t23_coinc_window + 1
        n2_incr = T_2 * T_3 * max(T_2, T_3) * slide_no
        n3_incr = T_2 * T_3 * (max(T_2, T_3) + T_23) * slide_no
        assert n2_incr < len(idx), "too many timeslides"
        idx[:, 1] = np.roll(idx[:, 1], n2_incr)
        idx[:, 2] = np.roll(idx[:, 2], n3_incr)
    return idx


@njit
def three_det_sum_idx_jit(t1, t2, t3, idx):
    """Most efficient if you are using numba and have the indexes."""
    temp = np.array([t1[i] + t2[j] + t3[k] for i, j, k in idx])
    return temp


@njit
def two_det_sum_idx_jit(t1, t2, idx):
    return np.array([t1[i] + t2[j] for i, j in idx])


def get_index_array_dtype(max_index):
    # reduce the size of the index array where possible
    bits_dtypes = [(8, np.uint8), (16, np.uint16),
                   (32, np.uint32), (64, np.uint64)]
    for bits, dtype in bits_dtypes:
        if 2**bits > max_index + 1:
            break
    return dtype


@njit
def get_indices_jit_3_ifo(tlen, t2_coinc_window,
                          t3_coinc_window, dtype=np.int64):
    idx = np.array(
        [
            [i, j, k]
            for i in range(tlen)
            for j in range(max(i - t2_coinc_window, 0), 
                           min(tlen, i + t2_coinc_window + 1))
            for k in range(max(i - t3_coinc_window, 0), 
                           min(tlen, i + t3_coinc_window + 1))
        ],
        dtype=dtype,
    )
    return idx


@njit
def get_indices_jit_2_ifo(tlen, t2_coinc_window, dtype=np.int64):
    idx = np.array(
        [
            [i, j]
            for i in range(tlen)
            for j in range(max(i - t2_coinc_window, 0), 
                           min(tlen, i + t2_coinc_window + 1))
        ],
        dtype=dtype,
    )
    return idx


def index_combinations(tlen, t2_coinc_window, t3_coinc_window, 
                       dtype=np.int64):
    """For three detectors, this is equivalent to calling 
    get_indices_jit_3_ifo, but is generally faster. For 
    two detectors this just calls get_indices_jit_2_ifo 
    directly.
    """
    if t3_coinc_window is None:
        return get_indices_jit_2_ifo(tlen, t2_coinc_window, dtype)
    elif 2 * max(t2_coinc_window, t3_coinc_window) > tlen:
        return get_indices_jit_3_ifo(tlen, t2_coinc_window, t3_coinc_window,
                                     dtype)
    # forgetting the tails at first
    largest_window = max(t2_coinc_window, t3_coinc_window)
    idx_1_mid = (
        np.arange(tlen - 2 * (t3_coinc_window + 1), dtype=dtype).repeat(
            (2 * t2_coinc_window + 1) * (2 * t3_coinc_window + 1)
        )
        + largest_window
        + 1
    )
    idx_2_mid = idx_1_mid + np.tile(
        np.arange(-t2_coinc_window, t2_coinc_window + 1, dtype=dtype).repeat(
            2 * t3_coinc_window + 1
        ),
        (tlen - 2 * (t3_coinc_window + 1)),
    )
    idx_3_mid = idx_1_mid + np.tile(
        np.arange(-t3_coinc_window, t3_coinc_window + 1, dtype=dtype),
        (tlen - 2 * (t3_coinc_window + 1)) * (2 * t2_coinc_window + 1),
    )
    idx_1_ends, idx_2_ends, idx_3_ends = get_indices_jit_3_ifo(
        2 * (t3_coinc_window + 1), t2_coinc_window, t3_coinc_window,
        dtype=dtype
    ).T
    # now get the tails
    n_start = int(len(idx_1_ends) / 2)
    idx_1 = np.concatenate(
        (
            idx_1_ends[:n_start],
            idx_1_mid,
            idx_1_ends[-n_start:] + tlen - (2 * t3_coinc_window + 2),
        )
    )
    idx_2 = np.concatenate(
        (
            idx_2_ends[:n_start],
            idx_2_mid,
            idx_2_ends[-n_start:] + tlen - (2 * t3_coinc_window + 2),
        )
    )
    idx_3 = np.concatenate(
        (
            idx_3_ends[:n_start],
            idx_3_mid,
            idx_3_ends[-n_start:] + tlen - (2 * t3_coinc_window + 2),
        )
    )
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
    snr_2_filt_rss = abs(snr_2_filt_rss) ** 2
    nifos = len(snr_2_filt_rss)
    dtype = get_index_array_dtype(np.max(idx))
    if nifos == 3:
        network_snr_sq = three_det_sum_idx_jit(
            snr_2_filt_rss[0], snr_2_filt_rss[1], snr_2_filt_rss[2],
            idx.astype(dtype)
        )
    elif nifos == 2:
        network_snr_sq = two_det_sum_idx_jit(
            snr_2_filt_rss[0], snr_2_filt_rss[1], idx.astype(dtype)
        )
    mask = network_snr_sq > threshold**2
    return np.sqrt(network_snr_sq[mask]), idx[mask]


def inner_complex(a, b):
    """Assumes a 2D array of detector SNRs where zeroth index selects a 
    detector.
    """
    return abs(np.sum(a * b.conjugate(), axis=0))


def snr_2_filter_and_threshold(snr_dom, snr_sub_perp, idx, threshold):
    """Calculate the 2 filter snr, and return values above threshold.

    Parameters
    ----------
    snr_dom : pycbc.TimeSeries
        _description_
    snr_sub_perp : pycbc.TimeSeries
        _description_
    idx : numpy.array
        _description_
    threshold : float
        SNR threshold

    Returns
    -------
    snr_2_filter: numpy.array
        _description_
    idx: numpy.array
    mask: boolean numpy.array
    """
    nifos = len(snr_dom)
    dom = np.array([snr_dom[i][idx[:, i]] for i in range(nifos)])
    sub = np.array([snr_sub_perp[i][idx[:, i]] for i in range(nifos)])
    network_rho_dom = np.linalg.norm(dom, axis=0)
    network_rho_sub_perp = inner_complex(dom, sub) / network_rho_dom
    snr_2_filter = np.sqrt(network_rho_sub_perp**2 + network_rho_dom**2)
    mask = snr_2_filter > threshold
    return snr_2_filter[mask], idx[mask], mask


@njit
def _check_time_idx_is_sorted_bool(time_idx):
    cond = [time_idx[i] <= time_idx[i + 1] for i in range(len(time_idx) - 1)]
    return cond


def check_time_idx_is_sorted(time_idx):
    """Wrapper for _check_time_idx_is_sorted_bool.

    Parameters
    ----------
    time_idx : numpy.array
        List of time indices.
    """
    cond = _check_time_idx_is_sorted_bool(time_idx)
    assert np.all(cond), "time_idx must be sorted"


def maximal_coinc_in_ifo(snrs, time_idx, check_sorted=True):
    """Choose the maximum network snr for each time point.
    Requires time_idx to be sorted.
    Parameters
    ----------
    snrs: numpy.array
        SNR-like values to maximize over.
    time_idx: numpy.array
        Sorted time indices.
    Returns
    -------
    i_max: numpy.array
        The indices that maximize the snr.
    """
    unsorted = False
    if check_sorted:
        try:
            check_time_idx_is_sorted(time_idx)
        except AssertionError:
            unsorted = True
            original_indices = np.arange(len(time_idx))
            sort_idx = np.argsort(time_idx)
            original_indices = original_indices[sort_idx]
            snrs = snrs[sort_idx]
            time_idx = time_idx[sort_idx]
    i_max = []
    j_start, j_end = 0, 0
    for c in np.unique(time_idx, return_counts=True)[1]:
        j_end += c
        i_max.append(np.argmax(snrs[j_start:j_end]) + j_start)
        j_start += c
    i_max = np.array(i_max)
    if unsorted:
        i_max = original_indices[i_max]
    return i_max


def maximize_snr_per_timepoint(snrs, det_idx, check_sorted=True):
    """Choose the maximum network snr for each time point in each detector.
    Assumes three detector network.
    Parameters
    ----------
    snrs: numpy.array
        SNR-like values to maximize over.
    det_idx: numpy.array2d
        Time indices. Must be sorted in zeroth detector.
    Returns
    -------
    snrs: numpy.array
        The maximum SNR-like values.
    det_idx: numpy.array2d
        Maximum SNR time indices. Sorted in zeroth detector time.
    i_max: numpy.array
        The indices that maximize the snr.
    """
    nifos = np.shape(det_idx)[-1]
    i_max = [maximal_coinc_in_ifo(snrs, det_idx[:, 0],
                                  check_sorted=check_sorted)]
    snrs = snrs[i_max[0]]
    det_idx = det_idx[i_max[0]]
    for i in range(1, nifos):
        sort_idx = np.argsort(det_idx[:, i])
        snrs = snrs[sort_idx]
        det_idx = det_idx[sort_idx]
        idx = maximal_coinc_in_ifo(snrs, det_idx[:, i], check_sorted=False)
        i_max.append(i_max[i - 1][sort_idx][idx])
        snrs = snrs[idx]
        det_idx = det_idx[idx]
    # sort so that ifo 0 times are again in order.
    sort_idx = np.argsort(det_idx[:, 0])
    return snrs[sort_idx], det_idx[sort_idx], i_max[-1][sort_idx]
