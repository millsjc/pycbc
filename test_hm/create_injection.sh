# to run this file: $ . inject.sh m1 m2 dist approx ra dec inc geocent_time(start-time) coa_phase psi flow
# m1 and m2 are in msuns and detector frame masses
# distance is in kpc (not Mpc!!)
# all angles are in degrees
# currently ignores spin in the injection

#source _inject.sh 150 30 100000 'IMRPhenomHMpseudoFourPN' 0 0 90 1187008882 0 0 20
#overhead Hanford:
source _inject.sh 150 30 100000 'IMRPhenomHMpseudoFourPN' 36.95072802 46.45514666 90 1187008882 0 0 20
#source _inject.sh 1.457423 1.299302 30000 'IMRPhenomD' 0 0 0 1187008882 0 0 20
