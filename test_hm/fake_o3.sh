pycbc_multi_inspiral_hm \
--fake-strain H1:aLIGOaLIGOO3LowT1800545 L1:aLIGOaLIGOO3LowT1800545 V1:aLIGOAdVO3LowT1800545 \
--injection-file H1:injection.xml L1:injection.xml V1:injection.xml \
--allow-zero-padding \
--psd-model H1:aLIGOaLIGOO3LowT1800545 L1:aLIGOaLIGOO3LowT1800545 V1:aLIGOAdVO3LowT1800545 \
--channel-name H1:DCH-CLEAN_STRAIN_C02 L1:DCH-CLEAN_STRAIN_C02 V1:Hrec_hoft_16384Hz  \
--gps-start-time 1187008870 --gps-end-time 1187008886 \
--strain-high-pass 15 --sample-rate 512 \
--segment-length 16 --segment-start-pad 4 \
--segment-end-pad 2 --allow-zero-padding  --taper-data 1 \
--low-frequency-cutoff 20 \
--timing-error 0.005 \
--approximant 'IMRPhenomHM:mtotal<4' 'IMRPhenomHM:else' \
--order -1 --snr-threshold 4. --cluster-window 0 \
--processing-scheme cpu \
--output TESTING-512.hdf \
--instruments H1 L1 V1 \
--bank-file GW190521_template.hdf \
--coinc-threshold 7.0 --verbose --num-timeslides 2

#--injection-file H1:injection.xml L1:injection.xml V1:injection.xml \
#--ra 0.644911865031637 \
#--dec 0.810795263791696 --trigger-time 1187008882 \
#--gps-start-time 1187008374 --gps-end-time 1187008886 \
#1187008860 --gps-end-time 1187008886
#H1:aLIGOaLIGOO3LowT1800545 L1:aLIGOaLIGOO3LowT1800545 V1:aLIGOAdVO3LowT1800545
