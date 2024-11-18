#!/bin/bash

# create global variables 
export DAK_PARAMS=$1
export DAK_RESULTS=$2
export USER_PARAMS=$3
export DAK_EVAL_NUM=$(
  grep "eval_id" $DAK_PARAMS | awk '{print $1}'
)

# directories for this evaluation
export WORKING_DIR=$(pwd)
export LAUNCH_DIR=$SLURM_SUBMIT_DIR
export TEMPLATE_DIR=""

# absolute path of the driver
driver_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# setup pre- and post-processing files
export PREPROCESS_FILE="${driver_dir}/pre_processor.sh"
export POSTPROCESS_FILE="${driver_dir}/post_processor.sh"

# setup size and executables for fluid solver
export M2C_SIZE=""
export M2C_EXE=""
export M2C_INPUT=""

# setup size and executables for solid solver
export AEROS_SIZE=""
export AEROS_EXE=""
export AEROS_INPUT=""

# Evaluation concurrency for dakota
export EVALUATION_CONCURRENCY=""

# get user variables
if [[ $USER_PARAMS == "" ]]
then
  echo "*** Error: Configuration file with paths to solver \
    executables, inputs, and resource requirements not provided.\
    Aborting ..."
  exit 1
fi
source $USER_PARAMS

# perform checks before proceeding
source $driver_dir/checks.sh

# --------------
# PRE-PROCESSING
# --------------

# copy all templates to the working directory
if [[ -e $TEMPLATE_DIR/fem.in.template ]]
then
  cp $TEMPLATE_DIR/fem.in.template \
    $WORKING_DIR/fem.in
else
  echo "*** Error: Could not find a template file for Aero-S input \
    file. Aborting ..."
  exit 1
fi
if [[ -e $TEMPLATE_DIR/input.st.template ]]
then
  cp $TEMPLATE_DIR/input.st.template \
    $WORKING_DIR/input.st
else
  echo "*** Error: Could not find a template file for M2C input \
    file. Aborting ..."
  exit 1
fi
if [[ -e $TEMPLATE_DIR/SphericalShock.txt.template ]]
then
  cp $TEMPLATE_DIR/SphericalShock.txt.template \
    $WORKING_DIR/SphericalShock.txt
else
  echo "*** Error: A template file for initial detonation profile not \
    provided. Aborting ..."
  exit 1
fi
if [[ -e $TEMPLATE_DIR/struct.geo.template ]]
then
  cp $TEMPLATE_DIR/struct.geo.template \
    $WORKING_DIR/struct.geo
else
  echo "*** Error: A template file for GMSH for creating mesh for each \
    design point was not provided. Aborting ..."
  exit 1
fi

# execute pre-processor in a sub-shell
if ! bash $PREPROCESS_FILE; then
  echo "*** Error: Failed at pre-processing stage for desing ${DAK_EVAL_NUM}. \
    Aborting ..."
  exit 1
fi

# ---------
# EXECUTION
# ---------

# find correct relative node position
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
  ($DAK_EVAL_NUM-1) % $EVALUATION_CONCURRENCY
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
  if ! bash $POSTPROCESS_FILE; then
    echo "*** Error: Failed at post-processing stage for desing \
      ${DAK_EVAL_NUM}. Aborting ..."
    exit 1
  fi

else
  # unsuccessfull evaluation
  printf "FAIL\n" > $DAK_RESULTS
fi
