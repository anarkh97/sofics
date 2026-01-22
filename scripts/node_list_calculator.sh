#!/bin/bash

# exit early when SLURM is not available and we are not a SLURM job
if ! command -v sbatch > /dev/null 2>&1 || [[ -z "${SLURM_JOB_ID:-}" ]]; then
  exit 0
fi

solid_size=$1
fluid_size=$2

#------------------------------------------------------------------------------
# Initialize quantities 
#------------------------------------------------------------------------------

# get all nodes assigned to dakota/ADOPT process
node_list=(`scontrol show hostnames`)

# total MPI tasks allocated for each node --- user's SLURM script
task_per_node=$(
  echo "$SLURM_TASKS_PER_NODE" |
  sed 's|^[^0-9]*\([0-9]\+\).*|\1|'
)

# total processors required for this process
total_proc=$((
  "$fluid_size"+"$solid_size"
))

# total number of compute nodes required by this process
applic_nodes=$((
  ("$total_proc"+"$task_per_node"-1) / 
  "$task_per_node"
))

# relative evaluation number (0, 1, ..., EVALUATION_CONCURRENCY-1)
relative_eval_num=$((
  ("$DAK_EVAL_NUM"-1) % "$EVALUATION_CONCURRENCY"
))

#------------------------------------------------------------------------------
# Find correct relative node position
#------------------------------------------------------------------------------

rnode_list=""

if [[ "$LOCAL_SCHEDULING" == 'DYNAMIC' ]]; then
  
  printf "*** Error: Dynamic local evalution tiling is not "
  printf "supported currently.\n"
  exit 1

  # TODO: Add this utility in future. Use files to block cores.

else

  # Static scheduling. Once an evaluation is complete its modulo is
  # launched. E.g., if evaluation 2 is completed and the concurrency 
  # is set to 10, evaluation 12 is launched next.

  start_node=$((
    ( "$relative_eval_num" * "$total_proc" ) / "$task_per_node"
  ))
  rnode_list="${node_list[$start_node]}"
  for i in $(seq 1 $((applic_nodes - 1)))
  do
    rnode_list="$host_list,${node_list[$((start_node+i))]}"
  done

fi

echo "$rnode_list"

# -----------------------------------------------------------------------------
# Exit assuming everything was successfull
# -----------------------------------------------------------------------------
exit 0
