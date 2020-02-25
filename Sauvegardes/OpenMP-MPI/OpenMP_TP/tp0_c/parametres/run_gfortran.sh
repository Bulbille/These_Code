# Execution monoprocesseur
./hello_world_mono.exe

# Execution parallele avec 1 thread
export OMP_NUM_THREADS=1
./hello_world_omp.exe

# Execution parallele avec 2 threads
export OMP_NUM_THREADS=2
./hello_world_omp.exe

# Execution parallele avec 4 threads
export OMP_NUM_THREADS=4
./hello_world_omp.exe

# Execution parallele avec 6 threads
export OMP_NUM_THREADS=6
./hello_world_omp.exe

# Execution parallele avec 8 threads
export OMP_NUM_THREADS=8
./hello_world_omp.exe
