#!/bin/bash

#SBATCH --job-name=adopt_true       # Job name
#SBATCH --nodes=5                   # Number of nodes
#SBATCH --ntasks-per-node=128       # Number of tasks per node
#SBATCH --time=144:00:00            # Time limit hrs:min:sec
#SBATCH --partition=normal_q        # Partition or queue name
#SBATCH --account=m2clab            # Cluster account
##SBATCH --exclude=tc[003,106,158-159,248] # Having problems w/ these nodes.
##SBATCH --exclusive=user

export UCX_LOG_LEVEL=error

# Change to the directory from which the job was submitted
cd $SLURM_SUBMIT_DIR

### Run
if [[ -e adopt.rst ]]
then
  ~/tinkercliffs/adopt/adopt -i adopt.in -o adopt.log -r adopt.rst
else
  ~/tinkercliffs/adopt/adopt -i adopt.in -o adopt.log -w adopt.rst
fi

exit;
