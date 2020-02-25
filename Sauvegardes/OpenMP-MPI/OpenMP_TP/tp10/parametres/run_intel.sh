# Execution monoprocesseur
./strassen_mono.exe

export OMP_STACKSIZE=512M

# Execution parallele avec 1 thread
export OMP_NUM_THREADS=1
./strassen_omp.exe

# Execution parallele avec 2 threads
export OMP_NUM_THREADS=2
./strassen_omp.exe

# Execution parallele avec 4 threads
export OMP_NUM_THREADS=4
./strassen_omp.exe

# Execution parallele avec 6 threads
export OMP_NUM_THREADS=6
./strassen_omp.exe

# Execution parallele avec 8 threads
export OMP_NUM_THREADS=8
./strassen_omp.exe
