#!/bin/bash

# create global variables 
export USER_PARAMS=$1
export DAK_PARAMS=$2
export DAK_RESULTS=$3
export DAK_EVAL_NUM=$(
  grep "eval_id" $DAK_PARAMS | awk '{print $1}'
)

# directories for this evaluation
export WORKING_DIR=$(pwd)
export LAUNCH_DIR=$SLURM_SUBMIT_DIR
export DRIVER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export TEMPLATE_DIR=""

# setup pre- and post-processing files
export PREPROCESS_FILE="${DRIVER_DIR}/pre_processor.sh"
export POSTPROCESS_FILE="${DRIVER_DIR}/post_processor.sh"

# executable for finite element mesher
export GMSH_EXE=""
export GMSH_INPUT=""

# setup size and executables for fluid solver
export M2C_SIZE=""
export M2C_EXE=""
export M2C_INPUT=""

# setup size and executables for solid solver
export AEROS_SIZE=""
export AEROS_EXE=""
export AEROS_INPUT=""

# detonation input for fluid solver
export SHOCK_INPUT=""

# Evaluation concurrency for dakota
export EVALUATION_CONCURRENCY=""

# get user variables
if [[ $USER_PARAMS == "" ]]
then
  printf "*** Error: Configuration file with paths to solver "
  printf "executables, inputs, and resource requirements not provided."
  printf "Aborting ...\n"
  exit 1
fi
source $USER_PARAMS

# perform checks before proceeding
source "${DRIVER_DIR}/checks.sh"

# --------------
# PRE-PROCESSING
# --------------

# copy all templates to the working directory
if [[ -e "$TEMPLATE_DIR/$AEROS_INPUT".template ]]; then
  cp "$TEMPLATE_DIR/$AEROS_INPUT".template \
    $WORKING_DIR/$AEROS_INPUT
else
  printf "*** Error: Could not find a template file for Aero-S input "
  printf "file. Aborting ...\n"
  exit 1
fi
if [[ -e "$TEMPLATE_DIR/$M2C_INPUT".template ]]; then
  cp "$TEMPLATE_DIR/$M2C_INPUT".template \
    $WORKING_DIR/$M2C_INPUT
else
  printf "*** Error: Could not find a template file for M2C input "
  printf "file. Aborting ...\n"
  exit 1
fi
if [[ -e "$TEMPLATE_DIR/$SHOCK_INPUT".template ]]; then
  cp "$TEMPLATE_DIR/$SHOCK_INPUT".template \
    $WORKING_DIR/$SHOCK_INPUT
else
  printf "*** Error: A template file for initial detonation profile not "
  printf "provided. Aborting ...\n"
  exit 1
fi
if [[ -e "$TEMPLATE_DIR/$GMSH_INPUT".template ]]; then
  cp "$TEMPLATE_DIR/$GMSH_INPUT".template \
    $WORKING_DIR/$GMSH_INPUT
else
  printf "*** Error: A template file for GMSH for creating mesh for each "
  printf "design point was not provided. Aborting ...\n"
  exit 1
fi

# execute pre-processor in a sub-shell
if ! bash $PREPROCESS_FILE; then
  printf "*** Error: Failed at pre-processing stage for design "
  printf "${DAK_EVAL_NUM}. Aborting ...\n"
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

printf "\033[34mLaunching Evaluation ${DAK_EVAL_NUM} on nodes "
printf "${host_list[*]}\033[0m\n"

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

# now we monitor the log file in a sub-shell
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
    printf "*** Error: Failed at post-processing stage for design "
    printf "${DAK_EVAL_NUM}. Aborting ...\n"
    exit 1
  fi

else
  # unsuccessfull evaluation
  printf "FAIL\n" > $DAK_RESULTS
fi
