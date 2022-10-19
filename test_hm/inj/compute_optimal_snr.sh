pycbc_optimal_snr \
--input-file bbh-3dets.hdf \
--output-file bbh-3dets_snr.hdf \
--psd-model H1:aLIGOaLIGOO3LowT1800545 L1:aLIGOaLIGOO3LowT1800545 V1:aLIGOAdVO3LowT1800545 \
--ifos "H1" "L1" "V1" \
--sample-rate 512 \
--seg-length 16 \
--f-low 20 \
--verbose


        # new_inj_table["network_optimal_snr"] = np.sqrt(np.sum(
        #     new_inj_table[f"{ifo}_optimal_snr"] for ifo in opts.ifos
        # ))