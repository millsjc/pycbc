INJFILE='mini-bbh-3dets-5mins.hdf'
#GPS_START=1267732868
#GPS_START=1267735956
GPS_START=1267736456
OUTFILE='MINI-TESTING-512-5mins-fakedata.hdf'
#BANKFILE='full-correct-tau-xhm-spinning-2filter-listharms.hdf'
#BANKFILE='splitbank/full-correct-tau-xhm-spinning-2filter-listharms_SPLITBANK_0.hdf'
BANKFILE='mini-bank.hdf'
h1chan=H1:GWOSC-4KHZ_R1_STRAIN
l1chan=L1:GWOSC-4KHZ_R1_STRAIN
v1chan=V1:GWOSC-4KHZ_R1_STRAIN
h1frame=H1_GWOSC_O3b_4KHZ_R1
l1frame=L1_GWOSC_O3b_4KHZ_R1
v1frame=V1_GWOSC_O3b_4KHZ_R1
WD=/work/cameron.mills/projects/pycbc/test_hm
#WD=/Users/camill/projects/pycbc/test_hm
pycbc_multi_inspiral_hm \
--fake-strain 'zeroNoise' \
--allow-zero-padding \
--psd-model H1:aLIGOaLIGOO3LowT1800545 L1:aLIGOaLIGOO3LowT1800545 V1:aLIGOAdVO3LowT1800545 \
--channel-name H1:${h1chan} L1:${l1chan} V1:${v1chan}  \
--gps-start-time $GPS_START --gps-end-time 1267736468 \
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
--coinc-threshold 7.0 --verbose --num-timeslides 0

#--injection-file $WD/inj/$INJFILE \
#--psd-estimation median --psd-segment-length 16 --psd-segment-stride 8 --psd-inverse-length 16 --psd-num-segments 63 \
#--frame-type H1:${h1frame} L1:${l1frame} V1:${v1frame} \

#--psd-model H1:aLIGOaLIGOO3LowT1800545 L1:aLIGOaLIGOO3LowT1800545 V1:aLIGOAdVO3LowT1800545 \
# --fake-strain 'zeroNoise' \
#--bank-file bank/full-correct-tau-xhm-spinning-2filter-listharms.hdf \
#--bank-file inj/bbh-3dets_bank.hdf \
#--bank-file GW190521_template.hdf \
#--injection-file H1:injection.xml L1:injection.xml V1:injection.xml \
#--ra 0.644911865031637 \
#--dec 0.810795263791696 --trigger-time 1187008882 \
#--gps-start-time 1187008374 --gps-end-time 1187008886 \
#1187008860 --gps-end-time 1187008886

