""" This module contains functions for calculating and manipulating higher harmonic
2-filter triggers.
"""
import numpy, numpy.random
import pycbc.waveform, pycbc.filter, pycbc.types, pycbc.psd, pycbc.fft, pycbc.conversions
import numpy as np


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

def orthogonalize_modes(h_modes, psd=None, f_lower=20., fh=2048.0, dominant_mode='22', subdom_modes=["33"]):
    """
    Orthogonalize a set of waveforms for a given PSD. 
    Returns waveforms orthogonal to the dominant mode.
    
    Parameters
    ----------
    h_modes: dict
        Dictionary of waveform modes
    psd: pycbc.types.FrequencySeries 
        psd to use when orthogonalizing
    f_lower: float
        Low frequency cutoff
    fh: float
        High frequency cutoff
    dominant_mode: str
        Loudest mode to orthogonalize w.r.t. e.g. "22"
    subdom_modes: list
        Other modes to be orthogonalized e.g. ["33", "44"]
    
    Returns
    -------
    h_perp: dict
        Waveform modes orthogonal to the dominant mode
    """
    if h_modes[dominant_mode].sigma > 0:
        h_dom = h_modes[dominant_mode] / h_modes[dominant_mode].sigma
        h_dom.sigma = h_modes[dominant_mode].sigma
        h_dom.norm = h_dom.zeta = 1
        h_perp = {dominant_mode: h_dom}
        for lm in subdom_modes:
            h_sub = h_modes[lm]
            if h_sub.sigma > 0:
                h_sub /= h_sub.sigma
                # generate the orthogonal waveform
                zeta = pycbc.filter.overlap_cplx(h_dom, h_sub, psd=psd, 
                    low_frequency_cutoff=f_lower, high_frequency_cutoff=fh, 
                    normalized=False) 
                norm = 1 / (numpy.sqrt(1 - numpy.abs(zeta) ** 2))
                h_p_lm = (h_sub - zeta * h_dom) * norm
                h_p_lm.zeta = zeta
                h_p_lm.norm = norm
                h_p_lm.sigma = h_sub.sigma * norm
                h_perp[lm] = h_p_lm
            else:
                print("No power in mode %s" % lm)
                h_perp[lm] = h
                h_perp[lm].sigma = 0
    else:
        print("No power in any mode, returning original waveform")
        h_perp = h_modes
    return h_perp