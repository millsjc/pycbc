DESC='MULTIINSPIRAL_HM_fakedata_512_512s_snr_10'
INJFILE="mini-bbh-3dets-512s_snr_10.hdf"
#GPS_START=1267732868
#GPS_START=1267735956
PAD=2
GPS_START=$((1267735956-${PAD}))
GPS_END=$((1267736468-${PAD}))
DURATION=$(($GPS_END-$GPS_START))
OUTFILE="H1L1V1-${DESC}-${GPS_START}-${DURATION}.hdf"
#BANKFILE='full-correct-tau-xhm-spinning-2filter-listharms.hdf'
#BANKFILE='splitbank/full-correct-tau-xhm-spinning-2filter-listharms_SPLITBANK_0.hdf'
BANKFILE='mini-bank.hdf'
h1chan=H1:GWOSC-4KHZ_R1_STRAIN
l1chan=L1:GWOSC-4KHZ_R1_STRAIN
v1chan=V1:GWOSC-4KHZ_R1_STRAIN
h1frame=H1_GWOSC_O3b_4KHZ_R1
l1frame=L1_GWOSC_O3b_4KHZ_R1
v1frame=V1_GWOSC_O3b_4KHZ_R1
#WD=/work/cameron.mills/projects/pycbc/test_hm
WD=/Users/camill/projects/pycbc/test_hm
pycbc_multi_inspiral_hm \
--fake-strain H1:aLIGOaLIGOO3LowT1800545 L1:aLIGOaLIGOO3LowT1800545 V1:aLIGOAdVO3LowT1800545 \
--allow-zero-padding \
--psd-model H1:aLIGOaLIGOO3LowT1800545 L1:aLIGOaLIGOO3LowT1800545 V1:aLIGOAdVO3LowT1800545 \
--channel-name H1:${h1chan} L1:${l1chan} V1:${v1chan}  \
--gps-start-time $GPS_START --gps-end-time $GPS_END \
--strain-high-pass 15 --sample-rate 512 \
--segment-length 512 --segment-start-pad 4 \
--segment-end-pad 2 --allow-zero-padding  --taper-data 1 \
--low-frequency-cutoff 20 \
--timing-error 0.005 \
--approximant 'IMRPhenomXHM:mtotal<4' 'IMRPhenomXHM:else' \
--order -1 --snr-threshold 4 --cluster-window 0 \
--processing-scheme cpu \
--output $WD/inspiral/$OUTFILE \
--bank-file $WD/bank/$BANKFILE \
--instruments H1 L1 V1 \
--coinc-threshold 7.0 --verbose \
--num-timeslides 0 \
--injection-file $WD/inj/$INJFILE
