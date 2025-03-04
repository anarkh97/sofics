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
export M2C_AUX=""

# setup size and executables for solid solver
export AEROS_SIZE=""
export AEROS_EXE=""
export AEROS_INPUT=""

# Evaluation concurrency for dakota
export EVALUATION_CONCURRENCY=""

# get user variables
if [[ $USER_PARAMS == "" ]]
then
  printf "*** Error: Configuration file with paths to solver "
  printf "executables, inputs, and resource requirements not provided."
  printf "Aborting ...\n"

  # NOTE: Here we do not write a results file and force dakota to fail.
  # This is a brute force approach for cases when user specifies a failure
  # capture method other than `abort` in dakota.

  exit 1
fi

source $USER_PARAMS

# perform checks before proceeding
if ! source "${DRIVER_DIR}/checks.sh"; then
  # NOTE: Here we do not write a results file and force dakota to fail.
  # This is a brute force approach for cases when user specifies a failure
  # capture method other than `abort` in dakota.

  exit 1
fi

# --------------
# PRE-PROCESSING
# --------------

# copy all templates to the working directory
template_error=0
if [[ -e "$TEMPLATE_DIR/${AEROS_INPUT}.template" ]]; then
  cp "$TEMPLATE_DIR/${AEROS_INPUT}.template" "$WORKING_DIR/$AEROS_INPUT"
else
  printf "*** Error: Could not find a template file for Aero-S input "
  printf "file (or template extension is missing).\n"
  template_error=$((template_error+1))
fi
if [[ -e "$TEMPLATE_DIR/${M2C_INPUT}.template" ]]; then
  cp "$TEMPLATE_DIR/${M2C_INPUT}.template" "$WORKING_DIR/$M2C_INPUT"
else
  printf "*** Error: Could not find a template file for M2C input "
  printf "file (or template extension is missing).\n"
  template_error=$((template_error+1))
fi

if [[ -e "$TEMPLATE_DIR/${GMSH_INPUT}.template" ]]; then
  cp "$TEMPLATE_DIR/${GMSH_INPUT}.template" "$WORKING_DIR/$GMSH_INPUT"
else
  printf "*** Error: Could not find a template file for Gmsh input "
  printf "file (or template extension is missing).\n"
  template_error=$((template_error+1))
fi

# copy any auxilary inputs that were provided for m2c
if [[ $M2C_AUX != "" ]]; then
  # get names of all auxilary inputs
  IFS=: read -r -a fluid_aux_inps <<< "$M2C_AUX"

  # remove empty fields (or spaces)
  fluid_aux_inps=("${fluid_aux_inps[@]// /}")

  for i in "${!fluid_aux_inps[@]}"; do
    if [[ -e "$TEMPLATE_DIR/${fluid_aux_inps[$i]}.template" ]]; then
      cp "$TEMPLATE_DIR/${fluid_aux_inps[$i]}.template" \
        "$WORKING_DIR/${fluid_aux_inps[$i]}"
    else
      printf "*** Error: Could not find a template file for auxilary "
      printf "input ${fluid_aux_inps[$i]} (or template extension is "
      printf "missing).\n"
      template_error=$((template_error+1))
    fi
  done
fi

if [[ $template_error -gt 0 ]]; then

  # Here we let dakota capture the failure and proceed 
  # based on user specification.

  printf "FAIL\n" > $DAK_RESULTS
  exit 0
fi

# execute pre-processor in a sub-shell
if ! bash $PREPROCESS_FILE; then
  printf "*** Error: Failed at pre-processing stage for design "
  printf "${DAK_EVAL_NUM}.\n"

  # Here we let dakota capture the failure and proceed 
  # based on user specification.

  printf "FAIL\n" > $DAK_RESULTS
  exit 0
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
printf "${host_list[*]} "

### run analysis in background
mpiexec --bind-to none \
	-n $M2C_SIZE \
        --host $host_list:$total_proc \
	$M2C_EXE $WORKING_DIR/$M2C_INPUT : \
        -n $AEROS_SIZE \
	--host $host_list:$total_proc \
        $AEROS_EXE $WORKING_DIR/$AEROS_INPUT 2>&1 \
	| tee $WORKING_DIR/log.out > /dev/null &

mpi_pid=$!
printf "with process id $mpi_pid.\033[0m\n"

# -------------------
# WAIT FOR COMPLETION
# ------------------

wait $mpi_pid

# ---------------
# POST-PROCESSING
# ---------------

# check which signal was recieved and handle errors here
# Alternatively, we could also query $? as it returns the
# exit status of mpiexec.
if grep -q "NORMAL TERMINATION" $WORKING_DIR/log.out
then
  # succesfull evaluation
  # execute post-processor in a sub-shell
  if ! bash $POSTPROCESS_FILE; then
    printf "*** Error: Failed at post-processing stage for design "
    printf "${DAK_EVAL_NUM}.\n"
    # Here we let dakota capture the failure and proceed 
    # based on user specification.
    
    printf "FAIL\n" > $DAK_RESULTS
    exit 0
  fi

else
  # unsuccessfull evaluation
  printf "FAIL\n" > $DAK_RESULTS
fi
