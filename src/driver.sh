#!/bin/bash

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
DAK_EVAL_NUM=$(
  grep "eval_id" $DAK_PARAMS | awk '{print $1}'
)

#------------------------------------------------------------------------------
# Internal variables (local to current shell)
#------------------------------------------------------------------------------
# Internal variables for ADOPT.
SOLVER_TYPE="TRUE"
NEIGHBORS=""

# Internal variables
ANALYSIS_SETUP_FILE="${DRIVER_DIR}/structure_setup.sh"
META_GENERATOR_FILE="${DRIVER_DIR}/meta_generator.sh"

# Internal drivers
TRUE_EVALUATION_DRIVER="${DRIVER_DIR}/fsi_driver.sh"
APPROX_EVALUATION_DRIVER="${DRIVER_DIR}/approx_fsi_driver.sh"
ERROR_EVALUATION_DRIVER="${DRIVER_DIR}/error_driver.sh"

#------------------------------------------------------------------------------
# Get user config in the CURRENT shell
#------------------------------------------------------------------------------
if [ -n $USER_CONFIG ]; then
  source $USER_CONFIG
fi

#------------------------------------------------------------------------------
# Process parameters and call the structural pre-processor
#------------------------------------------------------------------------------
source $ANALYSIS_SETUP_FILE

if [[ "$SOLVER_TYPE" == "APPROX" || "$SOLVER_TYPE" == "ERROR" ]]; then

  if ! source "${DRIVER_DIR}/fest_checks.sh"; then
    # NOTE: Here we do not write a results file and force dakota to fail.
    # This is a brute force approach for cases when user specifies a failure
    # capture method other than `abort` in dakota.
  
    exit 1
  fi

  if ! bash $METAGENERATOR_FILE $NEIGHBORS ; then
    printf "*** Error: Failed to generate metafile for "
    printf "ADOPT (design ${DAK_EVAL_NUM}).\n"

    printf "FAIL\n" > $DAK_RESULTS
    exit 0
  fi

fi

# ---------
# EXECUTION
# ---------

case "SOLVER_TYPE" in
  
  TRUE)
    bash $TRUE_EVALUATION_DRIVER ;;
  APPROX)
    bash $APPROX_EVALUATION_DRIVER ;;
  ERROR)
    bash $ERROR_EVALUATION_DRIVER ;;
  *)
    echo "*** Error: Unknown solver type \"${SOLVER_TYPE}\"."
    exit 1
    ;;

esac

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
