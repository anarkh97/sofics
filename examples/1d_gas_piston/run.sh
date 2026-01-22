#!/bin/bash

#SBATCH --job-name=dakota           # Job name
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=126       # Number of tasks per node
#SBATCH --time=02:00:00             # Time limit hrs:min:sec
#SBATCH --partition=normal_q        # Partition or queue name
#SBATCH --account=m2clab            # Cluster account
##SBATCH --exclusive=user

export UCX_LOG_LEVEL=error

# Change to the directory from which the job was submitted
cd $SLURM_SUBMIT_DIR

### Run
if [[ -e dakota.rst ]]
then
  dakota -i dakota.in -o dakota.log -r dakota.rst
else
  dakota -i dakota.in -o dakota.log -w dakota.rst
fi

exit;
