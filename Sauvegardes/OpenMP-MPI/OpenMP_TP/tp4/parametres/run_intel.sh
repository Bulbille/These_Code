# Execution monoprocesseur
./cg_mono.exe

# Choix du OMP_SCHEDULE
export OMP_SCHEDULE="STATIC"
export KMP_STACKSIZE=600m

# Execution parallele avec 1 thread
export OMP_NUM_THREADS=1
./cg_omp.exe

# Execution parallele avec 2 threads
export OMP_NUM_THREADS=2
./cg_omp.exe

# Execution parallele avec 4 threads
export OMP_NUM_THREADS=4
./cg_omp.exe

# Execution parallele avec 6 threads
export OMP_NUM_THREADS=6
./cg_omp.exe

# Execution parallele avec 8 threads
export OMP_NUM_THREADS=8
./cg_omp.exe
