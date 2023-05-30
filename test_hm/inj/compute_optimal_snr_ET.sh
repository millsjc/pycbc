FILENAME=$1
OUTFILE=`python -c 'print("{}_snr.hdf".format("'$FILENAME'".split(".")[0]))'`
# OUTFILE="./snrs.xml"
../../bin/pycbc_optimal_snr \
--input-file $FILENAME \
--output-file $OUTFILE \
--asd-file E1:/Users/camill/projects/eventhm/eventhm/data/asds/curves_May_2019/et_d.txt E2:/Users/camill/projects/eventhm/eventhm/data/asds/curves_May_2019/et_d.txt E3:/Users/camill/projects/eventhm/eventhm/data/asds/curves_May_2019/et_d.txt I1:/Users/camill/projects/eventhm/eventhm/data/asds/curves_May_2019/ce_3Hz_prepended.txt L1:/Users/camill/projects/eventhm/eventhm/data/asds/curves_May_2019/ce_3Hz_prepended.txt \
--ifos E1 E2 E3 I1 L1 \
--sample-rate 512 \
--seg-length 4 \
# --f-low 10 \
# --verbose

result=$(python <<EOF
import h5py, numpy as np
with h5py.File("$OUTFILE", "r+") as f:
    f["network_optimal_snr"] = np.sqrt(
        np.sum([np.array(f[f"optimal_snr_{ifo}"])**2
                for ifo in ["E1", "E2", "E3", "I1", "L1"]], 0))
EOF
)
