#!/bin/bash

# Add comment here about available directory variables

TEMPLATE_DIR=$LAUNCH_DIR/templates

EVALUATION_CONCURRENCY=$(
  grep "evaluation_concurrency" $LAUNCH_DIR/dakota.in |
  cut -d "=" -f 2 
)

# setup size and executables for fluid solver
M2C_SIZE=56
M2C_INPUT=input.st
M2C_AUX=SphericalShock.txt

# setup size and executables for solid solver
AEROS_SIZE=8
AEROS_INPUT=fem.in

# setup executable and input file for FE mesher
GMSH_EXE=gmsh
GMSH_INPUT=struct.geo
