#PBS -S /bin/ksh             # ligne obligatoire specifiant le shell: ksh ou csh
#PBS -q multi                # classe multiple
#PBS -N omp_tp6              # nom de la requete (defaut : nom du fichier)
#PBS -l cpunum_job=8         # reservation de 8 processeurs
#PBS -j o                    # concatene la sortie standard avec l erreur standard
#PBS -l cputim_prc=01:00:00  # temps max (HH:MM:SS)/ process 
#PBS -l cputim_job=01:05:00  # temps max (HH:MM:SS)/job
#PBS -l memsz_job=15gb        # memoire max /job 
 
# Echo des commandes passees 
set -x 

# On se met dans le repertoire de soumission
cd $PBS_O_WORKDIR

# Execution
./run.sh

