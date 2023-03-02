OUTDIR="../results/foundmissed/"
TRIGS="../inspiral/H1L1V1-MULTIINSPIRAL_HM_fakedata_512_512s_snr_10-1267735954-512.hdf"
INJFILE="../inj/mini-bbh-3dets-512s_snr_10.hdf"
/Users/camill/projects/pycbc/bin/hm/pycbc_hm_inj_finder -f $TRIGS -j $INJFILE -W 1 -c 'snr_2_filter' -o $OUTDIR -v
