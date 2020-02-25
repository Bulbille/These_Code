[ -z "$1" ] && echo "Usage: ./run_solutions.sh <exe>" && exit 1

export MKL_NUM_THREADS=1

# Execution monoprocesseur
timeout 60s ./$1_mono.exe

# Choix du OMP_SCHEDULE
export OMP_SCHEDULE="STATIC"

# Execution parallele avec 1 thread
export OMP_NUM_THREADS=1
timeout 60s ./$1_omp.exe

# Execution parallele avec 2 threads
export OMP_NUM_THREADS=2
timeout 60s ./$1_omp.exe

# Execution parallele avec 4 threads
export OMP_NUM_THREADS=4
timeout 60s ./$1_omp.exe

# Execution parallele avec 6 threads
export OMP_NUM_THREADS=6
timeout 60s ./$1_omp.exe

# Execution parallele avec 8 threads
export OMP_NUM_THREADS=8
timeout 60s ./$1_omp.exe
