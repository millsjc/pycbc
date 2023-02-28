#!/bin/bash

# delete previous submissions
rm submit.sh

# create modified versions of atlas_real_o3.sh for each injection set
for SNR in 8 9 10
do
sed "s/INJFILE='mini-bbh-3dets-512s.hdf'/INJFILE='mini-bbh-3dets-512s_snr_${SNR}.hdf'/g" atlas_real_o3.sh > atlas_real_o3_snr_${SNR}.sh
sed 's|real_o3|atlas_real_o3_snr_'${SNR}'.sh' ../submission.sub > atlas_real_o3_snr_${SNR}.sh
echo 'condor_submit ${PWD}/atlas_real_o3_snr_'${SNR}'.sub' >> submit.sh
done

# create a final modified version of atlas_real_o3.sh with no injections
sed "/INJFILE='mini-bbh-3dets-512s.hdf'/d; /--injection-file $WD/inj/$INJFILE/d" atlas_real_o3.sh > atlas_real_o3_noinj.sh
sed 's|real_o3|atlas_real_o3_noinj.sh' ../submission.sub > atlas_real_o3_noinj.sub
echo 'condor_submit ${PWD}/atlas_real_o3_noinj.sub' >> submit.sh

chmod +x atlas_real_o3_snr_*.sh
chmod +x atlas_real_o3_noinj.sh
