#!/bin/bash

# Add comment here about available directory variables

TEMPLATE_DIR=$LAUNCH_DIR/templatedir

EVALUATION_CONCURRENCY=$(
  grep "evaluation_concurrency" $LAUNCH_DIR/dakota.in |
  cut -d "=" -f 2 
)

# setup size and executables for fluid solver
M2C_SIZE=3
M2C_INPUT=input.st
M2C_AUX=StateCalculator.so

# setup size and executables for solid solver
AEROS_EXE="$LAUNCH_DIR/piston/piston"
AEROS_SIZE=1
AEROS_INPUT=spring.st
