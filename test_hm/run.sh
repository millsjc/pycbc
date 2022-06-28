#!/bin/bash
pycbc_multi_inspiral_hm \\
--verbose \\
--output test.hdf \\
--instruments 'L1,H1,V1' \\
--bank-file \\
--snr-threshold 4 \\
--newsnr-threshold 6 \\
--low-frequency-cutoff 20 \\
--approximant 'IMRPhenomXAS' \\
--cluster-method 'window' \\
--cluster-window 0.1 \\
--bank-veto-bank-file \\


