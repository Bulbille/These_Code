# Execution monoprocesseur
./poisson_mono.exe

# Choix du OMP_SCHEDULE
export OMP_SCHEDULE="STATIC"

# Execution parallele avec 1 thread
export OMP_NUM_THREADS=1
./poisson_omp.exe

# Execution parallele avec 2 threads
export OMP_NUM_THREADS=2
./poisson_omp.exe

# Execution parallele avec 4 threads
export OMP_NUM_THREADS=4
./poisson_omp.exe

# Execution parallele avec 6 threads
export OMP_NUM_THREADS=6
./poisson_omp.exe

# Execution parallele avec 8 threads
export OMP_NUM_THREADS=8
./poisson_omp.exe
