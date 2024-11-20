#!/bin/bash

# This utility simply performs checks for required variables
# for successfull evaluation of a design configuration.

PACKAGE_DIR=$(dirname $DRIVER_DIR)

if [[ $TEMPLATE_DIR == "" ]]; then
  printf "*** Warning: User did not provide a template directory. Using "
  printf "the default, which is the directory in which Dakota was executed.\n"
  TEMPLATE_DIR=$LAUNCH_DIR
fi

if [[ $GMSH_EXE == "" ]]; then
  printf "*** Error: Executable for GMSH (FE mesher) was not provided. "
  printf "Aborting ...\n"
  #printf "*** Warning: Executable for GMSH (FE mesher) was not provided. "
  #printf "Using default one ...\n"
  #if [[ ! -e $PACKAGE_DIR/install/bin/gmsh ]]; then
  #  printf "*** Error: No executable for solid solver was found. "
  #  printf "Aborting ...\n"
  #  exit 1
  #fi
  #GMSH_EXE=$PACKAGE_DIR/install/bin/gmsh
fi

if [[ $GMSH_INPUT == "" ]]; then
  printf "*** Error: GMSH (FE mesher) input file was not provided. "
  printf "Aborting ...\n"
  exit 1
fi

if [[ $M2C_SIZE == "" ]]; then
  printf "*** Error: Resource requirements for fluid solver were not "
  printf "specified. Aborting ...\n"
  exit 1
fi

if [[ $M2C_EXE == "" ]]; then
  printf "*** Warning: Executable for the fluid solver was not provided. "
  printf "Using default one ...\n"
  if [[ ! -e $PACKAGE_DIR/m2c/m2c ]]; then
    printf "*** Error: No executable for fluid solver was found. "
    printf "Aborting ...\n"
    exit 1
  fi
  M2C_EXE=$PACKAGE_DIR/m2c/m2c
fi

if [[ $M2C_INPUT == "" ]]; then
  printf "*** Error: Fluid input file was not provided. Aborting ...\n"
  exit 1
fi

if [[ $AEROS_SIZE == "" ]]; then
  printf "*** Error: Resource requirements for solid solver were not "
  printf "specified. Aborting ...\n"
  exit 1
fi

if [[ $AEROS_EXE == "" ]]; then
  printf "*** Warning: Executable for the solid solver was not provided. "
  printf "Using default one ...\n"
  if [[ ! -e $PACKAGE_DIR/aeros/bin/aeros ]]; then
    printf "*** Error: No executable for solid solver was found. "
    printf "Aborting ...\n"
    exit 1
  fi
  AEROS_EXE=$PACKAGE_DIR/aeros/bin/aeros
fi

if [[ $AEROS_INPUT == "" ]]; then
  printf "*** Error: Solid input file was not provided. Aborting ...\n"
  exit 1
fi

if [[ $SHOCK_INPUT == "" ]]; then
  printf "*** Error: File with initial detonation profile was not provided. "
  printf "Aborting ...\n"
  exit 1
fi

if [[ $EVALUATION_CONCURRENCY == "" ]]; then
  printf "*** Error: Evaluation concurrency needs to be defined in the "
  printf "configuration file as well. Should match dakota input file. "
  printf "Aborting ...\n"
  exit 1
fi
