../../bin/pycbc_optimal_snr \
--input-file bbh-3dets.hdf \
--output-file bbh-3dets_snr.hdf \
--psd-model H1:aLIGOaLIGOO3LowT1800545 L1:aLIGOaLIGOO3LowT1800545 V1:aLIGOAdVO3LowT1800545 \
--ifos "H1" "L1" "V1" \
--sample-rate 512 \
--seg-length 16 \
--f-low 20 \
--verbose

result=$(python <<EOF
import h5py, numpy as np
with h5py.File("bbh-3dets_snr.hdf", "r+") as f:
    f["network_optimal_snr"] = np.sqrt(
        np.sum([np.array(f[f"optimal_snr_{ifo}"])**2
                for ifo in ["H1", "L1", "V1"]], 0))
EOF
)
