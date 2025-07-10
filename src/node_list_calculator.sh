#!/bin/bash

rnode_list=$1
solid_size=$2
fluid_size=$3

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

# get all nodes assigned to dakota/ADOPT process
node_list=(`scontrol show hostnames`)

# calculate the relative node list for this process
host_list="${node_list[$relative_node]}"
for i in $(seq 1 $((applic_nodes - 1)))
do
  host_list="$host_list,
  ${node_list[$((relative_node+i))]}"
done
