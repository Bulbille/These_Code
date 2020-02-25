# Execution monoprocesseur
./fft_mono.exe

# Gestion de la taille de la memoire stack (valeur en KB)
export GOMP_STACKSIZE=10000000

# Choix du OMP_SCHEDULE
export OMP_SCHEDULE="STATIC"

# Execution parallele avec 1 thread
export OMP_NUM_THREADS=1
./fft_omp.exe

# Execution parallele avec 2 threads
export OMP_NUM_THREADS=2
./fft_omp.exe

# Execution parallele avec 4 threads
export OMP_NUM_THREADS=4
./fft_omp.exe

# Execution parallele avec 6 threads
export OMP_NUM_THREADS=6
./fft_omp.exe

# Execution parallele avec 8 threads
export OMP_NUM_THREADS=8
./fft_omp.exe
