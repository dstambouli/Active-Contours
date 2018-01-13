#!/bin/bash

# Slurm submission script,
# MPI job with Intel MPI/srun
# CRIHAN v 1.00 - Jan 2017
# support@criann.fr

# Not shared resources
#SBATCH --exclusive

# Partition (submission class)
#SBATCH --partition tcourt

# Job time (hh:mm:ss)
#SBATCH --time 00:20:00

# ----------------------------
# MPI tasks number
#SBATCH --ntasks 1

# MPI task maximum memory (MB)
#SBATCH --mem-per-cpu 3000
# ----------------------------

#SBATCH --mail-type ALL
# User e-mail address
##SBATCH --mail-user firstname.name@address.fr

# Compiler / MPI environments
# ------------------------------
module purge >& /dev/null
module load compilers/intel/2017
module load mpi/intelmpi/2017
# ------------------------------

LOGFILE=log.cv_tasks

# MPI code execution
srun cv | tee $LOGFILE
