#!/bin/bash

# This utility simply performs checks for required variables
# for successfull evaluation of a design configuration.

if [[ $M2C_SIZE == "" ]]; then
  echo "*** Error: Resource requirements for fluid solver were not specified. \
    Aborting ..."
  exit 1
fi

if [[ $M2C_EXE == "" ]]; then
  echo "*** Error: Executable for the fluid solver was not provided. Aborting ..."
  exit 1
fi

if [[ $M2C_INPUT == "" ]]; then
  echo "*** Error: Fluid input file was not provided. Aborting ..."
  exit 1
fi

if [[ $AEROS_SIZE == "" ]]; then
  echo "*** Error: Resource requirements for solid solver were not specified. Aborting ..."
  exit 1
fi

if [[ $AEROS_EXE == "" ]]; then
  echo "*** Error: Executable for the solid solver was not provided. Aborting ..."
  exit 1
fi

if [[ $AEROS_INPUT == "" ]]; then
  echo "*** Error: Solid input file was not provided. Aborting ..."
  exit 1
fi

if [[ $EVALUATION_CONCURRENCY == "" ]]; then
  echo "*** Error: Evaluation concurrency needs to be defined in the \
    configuration file as well. Should match dakota input file. \
    Aborting ..."
  exit 1
fi
