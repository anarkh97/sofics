#!/bin/bash

#SBATCH --job-name=dakota           # Job name
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=128       # Number of tasks per node
#SBATCH --time=00:10:00             # Time limit hrs:min:sec
#SBATCH --partition=normal_q        # Partition or queue name
#SBATCH --account=m2clab            # Cluster account

# Add any modules you might require. 
module load iimpi/2019b

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
