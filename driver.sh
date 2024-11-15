#!/bin/bash

# create global variables 
export DAK_PARAMS=$1
export DAK_RESULTS=$2
export DAK_EVAL_NUM=$(
  grep "eval_id" $DAK_PARAMS | awk '{print $1}'
)

# directories for this evaluation
export WORKING_DIR=$(pwd)
export LAUNCH_DIR=$SLURM_SUBMIT_DIR
export TEMPLATE_DIR="${LAUNCH_DIR}/templates"

# pre- and post-processing scripts
export PREPROCESS_FILE="pre_processor.sh"
export POSTPROCESS_FILE="post_processor.sh"

# setup size and executables for fluid solver
export M2C_SIZE=56
export M2C_EXE=~/tinkercliffs/m2c/m2c
export M2C_INPUT=$WORKING_DIR/input.st

# setup size and executables for solid solver
export AEROS_SIZE=8
export AEROS_EXE=~/tinkercliffs/FEMWorkingFoam/bin/aeros
export AEROS_INPUT=$WORKING_DIR/fem.in

# --------------
# PRE-PROCESSING
# --------------

# copy all templates to the working directory
cp $TEMPLATE_DIR/fem.in.template \
  $WORKING_DIR/fem.in
cp $TEMPLATE_DIR/input.st.template \
  $WORKING_DIR/input.st
cp $TEMPLATE_DIR/SphericalShock.txt.template \
  $WORKING_DIR/SphericalShock.txt
cp $TEMPLATE_DIR/struct.geo.template \
  $WORKING_DIR/struct.geo

# execute pre-processor in a sub-shell
bash $PREPROCESS_FILE

# ---------
# EXECUTION
# ---------

# find correct relative node position
eval_concurrency=$(
  grep "concurrency" $LAUNCH_DIR/dakota.in |
  cut -d "=" -f 2 
)
task_per_node=$(
  echo $SLURM_TASKS_PER_NODE |
  sed 's|^[^0-9]*\([0-9]\+\).*|\1|'
)
total_proc=$((
  $M2C_SIZE+$AEROS_SIZE
))
applic_nodes=$((
  ($total_proc+$task_per_node-1) / 
  $task_per_node
))
relative_eval_num=$((
  ($DAK_EVAL_NUM-1) % $eval_concurrency
))
relative_node=$((
  ( $relative_eval_num * $total_proc ) / $task_per_node
))

# get all nodes assigned to dakota process
node_list=(`scontrol show hostnames`)

# get host list for this process
host_list="${node_list[$relative_node]}"
for i in `seq 1 $(( $applic_nodes - 1))`
do
  host_list="$host_list,
  ${node_list[$((relative_node+i))]}"
done

printf '\033[34mLaunching Evaluation %s on nodes %s \033[0m\n'
  "${DAK_EVAL_NUM}" "${host_list[*]}"

### run analysis
mpiexec --bind-to none -n $M2C_SIZE \
  -host $host_list $M2C_EXE $M2C_INPUT : \
  -n $AEROS_SIZE -host $host_list \
  $AEROS_EXE $AEROS_INPUT &> $WORKING_DIR/log.out 

# -------------------
# WAIT FOR COMPLETION
# ------------------

# first we wait for log file to spawn
while [[ ! -e $WORKING_DIR/log.out ]]
do
  true
done

# now we monitor the log file
bash -c 'tail --pid=$$ -f "$1" |
  {
    sed -e "/Total Computation Time/q" \
        -e "/Found NAN in integrate/q" \
        -e "/ERROR/q" \
        -e "/Error/q" && kill $$;
  }' "MonitorShell" $WORKING_DIR/log.out

# ---------------
# POST-PROCESSING
# ---------------

# check which signal was recieved and handle errors here
if grep -q "Total Computation Time" $WORKING_DIR/log.out
then
  # succesfull evaluation
  # execute post-processor in a sub-shell
  bash $POSTPROCESS_FILE

else
  # unsuccessfull evaluation
  printf "FAIL\n" > $DAK_RESULTS
fi
