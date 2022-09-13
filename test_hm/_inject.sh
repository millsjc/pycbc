#!/bin/bash
# to run this file: $ . inject.sh m1 m2 dist approx ra dec inc geocent_time coa_phase psi flow
# m1 and m2 are in msuns and detector frame masses
# distance is in kpc (not Mpc!!)
# all angles are in degrees
# currently ignores spin in the injection
end_time=$(expr $8 + 10)

lalapps_inspinj --f-lower ${11} --disable-spin --waveform $4 --amp-order 0 --gps-start-time $8 --gps-end-time $end_time --t-distr uniform --l-distr fixed --longitude $5 --latitude $6 --m-distr fixMasses --fixed-mass1 $1 --fixed-mass2 $2 --i-distr fixed --fixed-inc $7 --coa-phase-distr fixed --fixed-coa-phase $9 --polarization ${10} --d-distr uniform --min-distance $3 --max-distance $3 --time-step 10 -o injection.xml
