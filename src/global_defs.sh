#!/bin/bash

#------------------------------------------------------------------------------
# Directory variables.
#------------------------------------------------------------------------------
export WORKING_DIR=$(pwd)
export LAUNCH_DIR=$SLURM_SUBMIT_DIR
export DRIVER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export TEMPLATE_DIR=""

#------------------------------------------------------------------------------
# Dakota Interface variables.
#------------------------------------------------------------------------------
export USER_CONFIG=""
export DAK_EVAL_NUM=""
export DAK_PARAMS=$DAKOTA_PARAMETERS_FILE
export DAK_RESULTS=$DAKOTA_RESULTS_FILE
export PREPROCESS_FILE="${DRIVER_DIR}/pre_processor.sh"
export POSTPROCESS_FILE="${DRIVER_DIR}/post_processor.sh"

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

# Evaluation concurrency for dakota
export EVALUATION_CONCURRENCY=""
