#                               -*- Mode: Sh -*- 
# job.sh ---  Job LoadLeveler de soumission de travaux sur IBM x3750
#
# Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
# Créé le         : Thu Jan 11 16:23:06 2001
# Dern. mod. par  : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
# Dern. mod. le   : Fri Aug 02 14:22:28 2013
#

# @ job_type = serial
# @ parallel_threads = 8
# @ wall_clock_limit = 300
# @ error  = omp_tp7.res
# @ output = omp_tp7.res
# @ class = cours
# @ queue

set -x

# Execution
./run.sh
