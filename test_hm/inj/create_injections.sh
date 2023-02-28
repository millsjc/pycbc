CONFIG=$1
SNR=(8 9 10)
OUTFILE='mini-bbh-3dets-512s.hdf'
PAD=2
GPS_START=$((1267735956-${PAD}))
GPS_END=$((1267736468-${PAD}))

pycbc_create_injections \
--gps-start-time ${GPS_START} \
--gps-end-time ${GPS_END} \
--time-step 10 --time-window 5 \
--seed 1234 \
--output-file $OUTFILE --force --verbose --config-files $CONFIG

#--gps-start-time 1267651175 \
#--gps-end-time 1267736468
# set in-plane spin to zero.
python remove_inplane_spin.py $OUTFILE

#generate identical bank
#python generate_identical_bank.py $OUTFILE

source compute_optimal_snr.sh $OUTFILE 

#if SNR variable is set then fix the SNR.
if [ -n "$SNR" ]; then
  if [ -n "$(declare -p SNR 2>/dev/null | grep 'declare -a')" ]; then
    # SNR is an array, fix each SNR value
    for snr_value in "${SNR[@]}"; do
      OUTFILE=`python -c 'print("{}_snr.hdf".format("'$FILENAME'".split(".")[0]))'`
      ./fix-snr.sh ${OUTFILE} ${snr_value}
    done
  else
    # SNR is a single value, fix it
    OUTFILE=`python -c 'print("{}_snr.hdf".format("'$FILENAME'".split(".")[0]))'`
    ./fix-snr.sh ${OUTFILE} $SNR
  fi
fi
