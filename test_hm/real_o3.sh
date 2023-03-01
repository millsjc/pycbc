DESC='MULTIINSPIRAL_HM_realdata'
INJFILE='mini-bbh-3dets-512s.hdf'
PAD=2
GPS_START=$((1267735956-${PAD}))
GPS_END=$((1267736468-${PAD}))
DURATION=$(($GPS_END-$GPS_START))
OUTFILE="H1L1V1-${DESC}-${GPS_START}-${DURATION}.hdf"
BANKFILE='mini-bank.hdf'
h1chan=H1:GWOSC-4KHZ_R1_STRAIN
l1chan=L1:GWOSC-4KHZ_R1_STRAIN
v1chan=V1:GWOSC-4KHZ_R1_STRAIN
h1frame=H1_GWOSC_O3b_4KHZ_R1
l1frame=L1_GWOSC_O3b_4KHZ_R1
v1frame=V1_GWOSC_O3b_4KHZ_R1
#FRAME_COMMAND="--frame-type H1:${h1frame} L1:${l1frame} V1:${v1frame}"
FRAME_COMMAND="--frame-files H1:./strain/H-${h1frame}-1267732480-4096.gwf L1:./strain/L-${l1frame}-1267732480-4096.gwf  V1:./strain/V-${v1frame}-1267732480-4096.gwf "
#WD=/work/cameron.mills/projects/pycbc/test_hm
WD=/Users/camill/projects/pycbc/test_hm
EXEC='pycbc_multi_inspiral_hm'


$EXEC \
--allow-zero-padding \
--psd-estimation median --psd-segment-length 16 --psd-segment-stride 8 --psd-inverse-length 16 --psd-num-segments 63 \
--channel-name H1:${h1chan} L1:${l1chan} V1:${v1chan}  \
$FRAME_COMMAND \
--gps-start-time $GPS_START --gps-end-time $GPS_END \
--pad-data $PAD \
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
