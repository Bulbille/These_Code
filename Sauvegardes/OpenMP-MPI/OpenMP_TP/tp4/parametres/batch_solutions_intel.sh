#                               -*- Mode: Sh -*- 
# job.sh ---  Job LoadLeveler de soumission de travaux sur IBM x3750
#
# Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
# Créé le         : Thu Jan 11 16:23:06 2001
# Dern. mod. par  : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
# Dern. mod. le   : Tue Jul 09 18:08:30 2013
#

# @ job_type = serial
# @ parallel_threads = 8
# @ wall_clock_limit = 1800
# @ error  = omp_tp4.res
# @ output = omp_tp4.res
# @ class = cours
# @ queue

# Echo des commandes passees
set -x

# On se met dans le repertoire de soumission
cd $LOADL_STEP_INITDIR

# Execution
for exec in $(ls *_omp.exe)
do
    exec=$(basename $exec _omp.exe)
    echo -n "Execution de $exec ... "
    ./run_solutions.sh $exec > omp_tp4_$exec.res 2>&1
    echo "done"
done
