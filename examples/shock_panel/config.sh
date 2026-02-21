#!/bin/bash

# Add comment here about available directory variables

TEMPLATE_DIR=$LAUNCH_DIR/templatedir

EVALUATION_CONCURRENCY=$(
  grep "evaluation_concurrency" $LAUNCH_DIR/adopt.in |
  cut -d "=" -f 2 
)

# setup size and executables for fluid solver
M2C_EXE=~/tinkercliffs/m2c_global/m2c
M2C_SIZE=63
M2C_INPUT=input.st

# setup size and executables for solid solver
AEROS_EXE=~/tinkercliffs/FEMWorkingFoam/bin/aeros
AEROS_SIZE=1
AEROS_INPUT=fem.in
