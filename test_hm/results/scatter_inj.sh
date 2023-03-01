for SNR in 8 9 10; do
    FOUNDMISSED_FILE="foundmissed/H1L1V1-MINI_TESTING_512_realdata_INJECTIONS_SNR_${SNR}_FOUNDMISSED-1267735954-512.h5"
    /Users/camill/projects/pycbc/bin/hm/pycbc_hm_scatter_inj_results --inj-file ${FOUNDMISSED_FILE} --bank-file ../bank/mini-bank.hdf --output-file ./figures/MINI-TESTING-512-realdata-INJECTIONS-missedfound_snr_${SNR}.png
done