#!/bin/bash

# Add comment here about available directory variables

TEMPLATE_DIR=$LAUNCH_DIR/templatedir

EVALUATION_CONCURRENCY=$(
  grep "evaluation_concurrency" $LAUNCH_DIR/adopt.in |
  cut -d "=" -f 2 
)

# setup size and executables for fluid solver
M2C_SIZE=30
M2C_INPUT=input.st

# setup size and executables for solid solver
AEROS_SIZE=1
AEROS_INPUT=fem.in
