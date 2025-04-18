#!/bin/bash

# Get the full-path to our driver script.
# NOTE: the "cd" and "pwd" commands are in a (...) 
#       i.e., they are executed in a sub-shell.
export DRIVER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#------------------------------------------------------------------------------
# Get all global variables in the CURRENT shell
#------------------------------------------------------------------------------
source "${DRIVER_DIR}/global_defs.sh"

#------------------------------------------------------------------------------
# Parse command-line in the SAME shell
#------------------------------------------------------------------------------
if ! source "${DRIVER_DIR}/command_line_parser.sh" "$@"
then
  exit 1
fi

#------------------------------------------------------------------------------
# Get current evaluation id
#------------------------------------------------------------------------------
set -o pipefail
DAK_EVAL_NUM=$(
  grep "eval_id" "$DAK_PARAMS" | awk '{print $1}'
)
if [ "$?" -ne 0 ]; then
  printf "*** Error: Could not get the evaluation ID from the "
  printf "parameters file. Required for SOFICS. Aborting ...\n"
  exit 1
fi
set +o pipefail

#------------------------------------------------------------------------------
# Internal variables (local to current shell)
#------------------------------------------------------------------------------
# Internal variables for ADOPT.
SOLVER_TYPE="TRUE"

# Internal scripts
ANALYSIS_SETUP_FILE="${DRIVER_DIR}/analysis_setup.sh"
META_GENERATOR_FILE="${DRIVER_DIR}/meta_generator.sh"

# Internal drivers
TRUE_EVALUATION_DRIVER="${DRIVER_DIR}/fsi_driver.sh"
APPROX_EVALUATION_DRIVER="${DRIVER_DIR}/approx_fsi_driver.sh"
ERROR_EVALUATION_DRIVER="${DRIVER_DIR}/error_driver.sh"

#------------------------------------------------------------------------------
# Get user config in the CURRENT shell
#------------------------------------------------------------------------------
if [ -n "$USER_CONFIG" ]; then
  source "$USER_CONFIG"
fi

#------------------------------------------------------------------------------
# Process parameters and call the structural pre-processor
#------------------------------------------------------------------------------
source "$ANALYSIS_SETUP_FILE"

#------------------------------------------------------------------------------
# Execute analysis in a sub-shell.
#------------------------------------------------------------------------------
case "$SOLVER_TYPE" in
  
  TRUE)
    bash "$TRUE_EVALUATION_DRIVER" ;;
  APPROX)
    bash "$APPROX_EVALUATION_DRIVER" ;;
  ERROR)
    bash "$ERROR_EVALUATION_DRIVER" ;;

esac

#------------------------------------------------------------------------------
# Post-processing
#------------------------------------------------------------------------------

# check which signal was recieved and handle errors here
# Alternatively, we could also query $? as it returns the
# exit status of mpiexec.
if grep -q "NORMAL TERMINATION" "$WORKING_DIR/log.out"
then
  # succesfull evaluation
  if [[ "$SOLVER_TYPE" != "ERROR" ]]; then
    
    # execute post-processor in a sub-shell
    if ! bash "$POSTPROCESS_FILE"; then
      printf "*** Error: Failed at post-processing stage for design "
      printf "%s.\n" "${DAK_EVAL_NUM}"
      # Here we let dakota capture the failure and proceed 
      # based on user specification.
      
      printf "FAIL\n" > "$WORKING_DIR/$DAK_RESULTS"
      exit 0
    fi

  else

    NRMSE=$(
      grep "Mean Squared Error" "$WORKING_DIR/log.out" |
      sed 's/.*: //; s/\.$//'
    )

    printf "    %s  MSE\n" "$NRMSE" >> "$WORKING_DIR/$DAK_RESULTS"
  fi
else
  # unsuccessfull evaluation
  printf "FAIL\n" > "$WORKING_DIR/$DAK_RESULTS"
fi

#------------------------------------------------------------------------------
# Exit analysis driver.
#------------------------------------------------------------------------------
exit 0
