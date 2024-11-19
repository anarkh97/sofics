#!/bin/bash

# This utility simply performs checks for required variables
# for successfull evaluation of a design configuration.

if [[ $GMSH_EXE == "" ]]; then
  printf "*** Error: Executable for GMSH (FE mesher) was not provided. "
  printf "Aborting ...\n"
  exit 1
fi

if [[ $M2C_SIZE == "" ]]; then
  printf "*** Error: Resource requirements for fluid solver were not "
  printf "specified. Aborting ...\n"
  exit 1
fi

if [[ $M2C_EXE == "" ]]; then
  printf "*** Error: Executable for the fluid solver was not provided. "
  printf "Aborting ...\n"
  exit 1
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
  printf "*** Error: Executable for the solid solver was not provided. "
  printf "Aborting ...\n"
  exit 1
fi

if [[ $AEROS_INPUT == "" ]]; then
  printf "*** Error: Solid input file was not provided. Aborting ...\n"
  exit 1
fi

if [[ $EVALUATION_CONCURRENCY == "" ]]; then
  printf "*** Error: Evaluation concurrency needs to be defined in the "
  printf "configuration file as well. Should match dakota input file. "
  printf "Aborting ...\n"
  exit 1
fi
