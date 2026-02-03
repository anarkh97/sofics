#!/bin/bash

# TODO: These could be converted to associative arrays.
#       Such data structures are useful for creating 
#       arrays of type, array[key] = value.
#       However, there is no easy way to export arrays.

#------------------------------------------------------------------------------
# Directory variables.
#------------------------------------------------------------------------------
export WORKING_DIR=$(pwd)
# when SLURM is not available and we are not a SLURM job
if ! command -v sbatch > /dev/null 2>&1 || [[ -z "${SLURM_JOB_ID:-}" ]]; then
  export LAUNCH_DIR=$(dirname "$WORKING_DIR")
else
  export LAUNCH_DIR=$SLURM_SUBMIT_DIR
fi
export TEMPLATE_DIR=""
export RESULTS_DIR=""

#------------------------------------------------------------------------------
# Dakota Interface variables.
#------------------------------------------------------------------------------
export USER_CONFIG=""
export DAK_EVAL_NUM=""
export DAK_PARAMS=$DAKOTA_PARAMETERS_FILE
export DAK_RESULTS=$DAKOTA_RESULTS_FILE
export PREPROCESS_FILE="${SOFICS_BIN}/pre_processor.sh"
export POSTPROCESS_FILE="${SOFICS_BIN}/post_processor.sh"

#------------------------------------------------------------------------------
# User config variables.
#------------------------------------------------------------------------------
# executable for finite element mesher
export GMSH_EXE=""
export GMSH_INPUT=""

# setup size and executables for fluid solver
export M2C_SIZE=""
export M2C_EXE=""
export M2C_INPUT=""
export M2C_AUX=""

# setup size and executables for approximate fluid solver
export FEST_SIZE=""
export FEST_EXE=""
export FEST_INPUT=""
export FEST_ERROR_INPUT=""

# setup size and executables for solid solver
export AEROS_SIZE=""
export AEROS_EXE=""
export AEROS_INPUT=""

# local scheduling option
export LOCAL_SCHEDULING=""

# evaluation concurrency for dakota
export EVALUATION_CONCURRENCY=""

#------------------------------------------------------------------------------
# ADOPT specific variables.
#------------------------------------------------------------------------------
export META_FILE=""
export META_SURFACE_FILE=""
export META_SOLUTION_FILE=""

