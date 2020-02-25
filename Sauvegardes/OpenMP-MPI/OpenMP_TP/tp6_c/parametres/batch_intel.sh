#                               -*- Mode: Sh -*- 
# job.sh ---  Job LoadLeveler de soumission de travaux sur IBM x3750
#
# Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
# Créé le         : Thu Jan 11 16:23:06 2001
# Dern. mod. par  : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
# Dern. mod. le   : Thu Aug 01 10:02:04 2013
#

# @ job_type = serial
# @ parallel_threads = 32
# @ wall_clock_limit = 300
# @ error  = omp_tp6.res
# @ output = omp_tp6.res
# @ class = cours
# @ queue

# Echo des commandes passees
set -x

# On se met dans le repertoire de soumission
cd $LOADL_STEP_INITDIR

# Execution
./run.sh
