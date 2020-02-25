#                               -*- Mode: Sh -*- 
# job.sh ---  Job LoadLeveler de soumission de travaux sur IBM SP
#
# Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
# Créé le         : Thu Jan 11 16:23:06 2001
# Dern. mod. par  : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
# Dern. mod. le   : Wed Jun 12 11:20:14 2013
#

# @ job_type = serial
# @ resources = ConsumableCpus(8)
# @ job_cpu_limit = 300
# @ error  = bi_cgstab.res
# @ output = bicgstab.res
# @ shell = /bin/ksh
# @ class = cours
# @ queue

set -x

# Compilation version monoprocesseur
make clean
make mono
# Execution monoprocesseur
./prod_mat

# Compilation version parallelisee avec OpenMP
make clean
make para
# Execution parallele avec 2 threads
export OMP_NUM_THREADS=2
export OMP_SCHEDULE="STATIC"
./prod_mat
# Execution parallele avec 4 threads
export OMP_NUM_THREADS=4
export OMP_SCHEDULE="STATIC"
./prod_mat
# Execution parallele avec 6 threads
export OMP_NUM_THREADS=6
export OMP_SCHEDULE="STATIC"
./prod_mat
# Execution parallele avec 8 threads
export OMP_NUM_THREADS=8
export OMP_SCHEDULE="STATIC"
./prod_mat

make clean
