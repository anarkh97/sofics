#!/bin/bash

# Add comment here about available directory variables

TEMPLATE_DIR=$LAUNCH_DIR/templates

# setup size and executables for fluid solver
M2C_SIZE=56
M2C_EXE=~/tinkercliffs/m2c/m2c
M2C_INPUT=$WORKING_DIR/input.st

# setup size and executables for solid solver
AEROS_SIZE=8
AEROS_EXE=~/tinkercliffs/FEMWorkingFoam/bin/aeros
AEROS_INPUT=$WORKING_DIR/fem.in
