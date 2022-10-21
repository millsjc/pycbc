OUTFILE='bbh-3dets-1hour.hdf'
GPS_START=1267732868
pycbc_create_injections \
--gps-start-time ${GPS_START} \
--gps-end-time 1267736468 \
--time-step 16 --time-window 5 \
--seed 1234 \
--output-file $OUTFILE --force --verbose --config-files injections.ini

#--gps-start-time 1267651175 \
#--gps-end-time 1267736468
# set in-plane spin to zero.
python remove_inplane_spin.py $OUTFILE

#generate identical bank
python generate_identical_bank.py $OUTFILE

source compute_optimal_snr.sh  
