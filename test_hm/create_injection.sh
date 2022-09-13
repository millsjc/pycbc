# to run this file: $ . inject.sh m1 m2 dist approx ra dec inc geocent_time(start-time) coa_phase psi flow
# m1 and m2 are in msuns and detector frame masses
# distance is in kpc (not Mpc!!)
# all angles are in degrees
# currently ignores spin in the injection

source _inject.sh 150 90 3000000 'IMRPhenomDpseudoFourPN' 0 0 60 1187008882 0 0 20
#source _inject.sh 1.457423 1.299302 30000 'IMRPhenomD' 0 0 0 1187008882 0 0 20
