INJFILE='bbh-3dets-1hour.hdf'
GPS_START=1267732868
OUTFILE='TESTING-512-1hour-zeronoise.hdf'
BANKFILE='full-correct-tau-xhm-spinning-2filter-listharms.hdf'
#WD=/work/cameron.mills/projects/pycbc/test_hm
#WD=/Users/camill/projects/pycbc/test_hm
pycbc_multi_inspiral_hm \
--fake-strain 'zeroNoise' \
--injection-file $WD/inj/$INJFILE \
--allow-zero-padding \
--psd-model H1:aLIGOaLIGOO3LowT1800545 L1:aLIGOaLIGOO3LowT1800545 V1:aLIGOAdVO3LowT1800545 \
--channel-name H1:DCH-CLEAN_STRAIN_C02 L1:DCH-CLEAN_STRAIN_C02 V1:Hrec_hoft_16384Hz  \
--gps-start-time $GPS_START --gps-end-time 1267736468 \
--strain-high-pass 15 --sample-rate 512 \
--segment-length 512 --segment-start-pad 4 \
--segment-end-pad 2 --allow-zero-padding  --taper-data 1 \
--low-frequency-cutoff 20 \
--timing-error 0.005 \
--approximant 'IMRPhenomXHM:mtotal<4' 'IMRPhenomXHM:else' \
--order -1 --snr-threshold 0.01 --cluster-window 0 \
--processing-scheme cpu \
--output $WD/$OUTFILE \
--bank-file $WD/bank/$BANKFILE \
--instruments H1 L1 V1 \
--coinc-threshold 7.0 --verbose --num-timeslides 0

#--bank-file bank/full-correct-tau-xhm-spinning-2filter-listharms.hdf \
#--bank-file inj/bbh-3dets_bank.hdf \
#--bank-file GW190521_template.hdf \
#--injection-file H1:injection.xml L1:injection.xml V1:injection.xml \
#--ra 0.644911865031637 \
#--dec 0.810795263791696 --trigger-time 1187008882 \
#--gps-start-time 1187008374 --gps-end-time 1187008886 \
#1187008860 --gps-end-time 1187008886
#H1:aLIGOaLIGOO3LowT1800545 L1:aLIGOaLIGOO3LowT1800545 V1:aLIGOAdVO3LowT1800545
