TRIG_FILE="../inspiral/H1L1V1-MINI_TESTING_512_realdata_NOINJECTIONS_3timeslides-1267735954-512.hdf"
OUTDIR="./"
pycbc_grb_trig_cluster --time-window 1 --rank-column "snr_2_filter" --trig-file $TRIG_FILE --output-dir $OUTDIR --verbose