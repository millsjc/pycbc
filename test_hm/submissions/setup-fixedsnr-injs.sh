#!/bin/bash

# delete previous submissions
rm submit.sh

# create modified versions of atlas_real_o3.sh for each injection set
for SNR in 8 9 10
do
sed "s/INJFILE='mini-bbh-3dets-512s.hdf'/INJFILE='mini-bbh-3dets-512s_snr_${SNR}.hdf'/g" scripts/atlas_real_o3.sh > scripts/atlas_real_o3_snr_${SNR}.sh
sed -i.bak "s#OUTFILE='MINI-TESTING-512-realdata-INJECTIONS-timeslides.hdf'#OUTFILE='MINI-TESTING-512-realdata-INJECTIONS-SNR-${SNR}.hdf'#g" scripts/atlas_real_o3_snr_${SNR}.sh
sed 's|real_o3|atlas_real_o3_snr_'${SNR}'|g' ../submission.sub > scripts/atlas_real_o3_snr_${SNR}.sub
echo 'condor_submit ${PWD}/scripts/atlas_real_o3_snr_'${SNR}'.sub' >> submit.sh
done

NUM_SLIDES=3
# create a final modified version of atlas_real_o3.sh with no injections and with timeslides
sed "s|INJFILE='mini-bbh-3dets-512s.hdf'||g" scripts/atlas_real_o3.sh > scripts/atlas_real_o3_noinj.sh
sed -i.bak 's|--injection-file $WD/inj/$INJFILE||g' scripts/atlas_real_o3_noinj.sh
sed -i.bak 's|--num-timeslides 0|--num-timeslides 3|g' scripts/atlas_real_o3_noinj.sh
sed -i.bak "s#OUTFILE='MINI-TESTING-512-realdata-INJECTIONS-timeslides.hdf'#OUTFILE='MINI-TESTING-512-realdata-NOINJECTIONS-${NUM_SLIDES}timeslides.hdf'#g" scripts/atlas_real_o3_noinj.sh
sed 's|real_o3|atlas_real_o3_noinj|g' ../submission.sub > scripts/atlas_real_o3_noinj.sub
echo 'condor_submit ${PWD}/scripts/atlas_real_o3_noinj.sub' >> submit.sh

chmod +x scripts/atlas_real_o3_snr_*.sh
chmod +x scripts/atlas_real_o3_noinj.sh
