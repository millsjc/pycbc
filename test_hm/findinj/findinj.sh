for SNR in 8 9 10; do
    TRIGS="../inspiral/H1L1V1-MINI_TESTING_512_realdata_INJECTIONS_SNR_${SNR}-1267735954-512.hdf"
    INJFILE="../inj/fixed-snr-mini-bbh-3dets-512s_snr_${SNR}.hdf"
    /Users/camill/projects/pycbc/bin/hm/pycbc_hm_inj_finder -f $TRIGS -j $INJFILE -W 1 -c 'snr_2_filter' -v
done
