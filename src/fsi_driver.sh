#!/bin/bash

#------------------------------------------------------------------------------
# Calculate node (list) for the current evaluation.
#------------------------------------------------------------------------------

# find correct relative node position
task_per_node=$(
  echo "$SLURM_TASKS_PER_NODE" |
  sed 's|^[^0-9]*\([0-9]\+\).*|\1|'
)
total_proc=$((
  "$M2C_SIZE"+"$AEROS_SIZE"
))
applic_nodes=$((
  ("$total_proc"+"$task_per_node"-1) / 
  "$task_per_node"
))
relative_eval_num=$((
  ("$DAK_EVAL_NUM"-1) % "$EVALUATION_CONCURRENCY"
))
relative_node=$((
  ( "$relative_eval_num" * "$total_proc" ) / "$task_per_node"
))

# get all nodes assigned to dakota process
node_list=(`scontrol show hostnames`)

# get host list for this process
host_list="${node_list[$relative_node]}"
for i in $(seq 1 $((applic_nodes - 1)))
do
  host_list="$host_list,
  ${node_list[$((relative_node+i))]}"
done

### run analysis in background
mpiexec --bind-to none \
	-n "$M2C_SIZE" \
  --host "$host_list":"$total_proc" \
	"$M2C_EXE" "$WORKING_DIR/$M2C_INPUT" : \
  -n "$AEROS_SIZE" \
	--host "$host_list":"$total_proc" \
  "$AEROS_EXE" "$WORKING_DIR/$AEROS_INPUT" 2>&1 \
	| tee "$WORKING_DIR/log.out" > /dev/null &

mpi_pid=$!

printf "\033[34mLaunching Evaluation %s on nodes " "${DAK_EVAL_NUM}"
printf "%s with process id %s.\033[0m\n" "${host_list[*]}" "$mpi_pid"

# -----------------------------------------------------------------------------
# Wait for completion
# -----------------------------------------------------------------------------

wait "$mpi_pid"

# -----------------------------------------------------------------------------
# Exit assuming everything was successfull
# -----------------------------------------------------------------------------
exit 0
