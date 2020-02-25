# @ job_type = serial
# @ parallel_threads = 8
# @ wall_clock_limit = 300
# @ error  = omp_tp0.res
# @ output = omp_tp0.res
# @ class = cours
# @ queue

# Echo des commandes passees
set -x

# On se met dans le repertoire de soumission
cd $LOADL_STEP_INITDIR

# Execution
./run.sh
