#!/bin/bash

# This utility simply performs checks for required variables
# for successfull evaluation of a design configuration.

BUILD_DIR=$(dirname "$DRIVER_DIR")
PACKAGE_DIR=$BUILD_DIR/packages

if [ -z "$TEMPLATE_DIR" ]; then
  printf "*** Warning: User did not provide a template "
  printf "directory. Using the default, which is the "
  printf "directory in which Dakota was executed.\n"
  TEMPLATE_DIR=$LAUNCH_DIR
fi

#------------------------------------------------------------------------------
# Executables checks. (Assign defaults whenever possible)
#------------------------------------------------------------------------------

# AN: Handled later.
#if [ -z "$GMSH_EXE" ]; then
#  printf "*** Error: Executable for GMSH (FE mesher) was not provided. "
#  printf "Aborting ...\n"
#  exit 1
#  #printf "*** Warning: Executable for GMSH (FE mesher) was not provided. "
#  #printf "Using default one ...\n"
#  #if [[ ! -e $PACKAGE_DIR/install/bin/gmsh ]]; then
#  #  printf "*** Error: No executable for solid solver was found. "
#  #  printf "Aborting ...\n"
#  #  exit 1
#  #fi
#  #GMSH_EXE=$PACKAGE_DIR/install/bin/gmsh
#elif [ -n "$GMSH_EXE" ] && [[ ! -x "$GMSH_EXE" ]]; then
#  printf "*** Error: \"${GMSH_EXE}\" is not a valid executable. "
#  printf "Aborting ...\n"
#  exit 1
#fi

if [ -z "$M2C_EXE" ]; then
  printf "*** Warning: Executable for M2C was not provided. "
  printf "Using default one ...\n"
  if [[ ! -x "$PACKAGE_DIR/m2c/m2c" ]]; then
    printf "*** Error: Could not find a valid executable for M2C. "
    printf "Aborting ...\n"
    exit 1
  fi
  M2C_EXE="$PACKAGE_DIR/m2c/m2c"
elif [ -n "$M2C_EXE" ] && [[ ! -x "$M2C_EXE" ]]; then
  printf "*** Error: \"%s\" is not a valid executable. " "$M2C_EXE"
  printf "Aborting ...\n"
  exit 1
fi

if [ -z "$AEROS_EXE" ]; then
  printf "*** Warning: Executable for AERO-S was not provided. "
  printf "Using default one ...\n"
  if [[ ! -x "$PACKAGE_DIR/aeros/bin/aeros" ]]; then
    printf "*** Error: Could not find a valid executable for AERO-S. "
    printf "Aborting ...\n"
    exit 1
  fi
  AEROS_EXE="$PACKAGE_DIR/aeros/bin/aeros"
elif [ -n "$AEROS_EXE" ] && [[ ! -x "$AEROS_EXE" ]]; then
  printf "*** Error: \"%s\" is not a valid executable. " "${AEROS_EXE}"
  printf "Aborting ...\n"
  exit 1
fi

#------------------------------------------------------------------------------
# Size checks. (Assign defaults whenever possible)
#------------------------------------------------------------------------------

if [ -z "$M2C_SIZE" ]; then
  printf "*** Error: Resource requirements for M2C were not "
  printf "specified. Aborting ...\n"
  exit 1
fi

if [ -z "$AEROS_SIZE" ]; then
  printf "*** Error: Resource requirements for AERO-S were not "
  printf "specified. Aborting ...\n"
  exit 1
fi

#------------------------------------------------------------------------------
# Input file checks. (Assign defaults whenever possible)
#------------------------------------------------------------------------------

template_error=0
#if [[ ! -f "$TEMPLATE_DIR/${GMSH_INPUT}.template" ]]; then
#  printf "*** Error: Could not find a template file for Gmsh input "
#  printf "file (or template extension is missing).\n"
#  template_error=$((template_error+1))
#fi

if [[ ! -f "$TEMPLATE_DIR/${AEROS_INPUT}.template" ]]; then
  printf "*** Error: Could not find the template file for AEROS input "
  printf "file (or template extension is missing). Aborting ...\n"
  template_error=$((template_error+1))
fi

if [[ ! -f "$TEMPLATE_DIR/${M2C_INPUT}.template" ]]; then
  printf "*** Error: Could not find a template file for M2C input "
  printf "file (or template extension is missing).\n"
  template_error=$((template_error+1))
fi

# auxilary inputs for m2c
if [ -n "$M2C_AUX" ]; then
  # convert ":" separated string into array
  IFS=: read -r -a fluid_aux_inps <<< "$M2C_AUX"
  for i in "${!fluid_aux_inps[@]}"; do
    if [[ ! -f "$TEMPLATE_DIR/${fluid_aux_inps[$i]}.template" ]]; then
      printf "*** Error: Could not find a template file for auxilary "
      printf "input %s (or template extension is " "${fluid_aux_inps[$i]}"
      printf "missing).\n"
      template_error=$((template_error+1))
    fi
  done
fi

# let driver handle these errors.
if [[ $template_error -gt 0 ]]; then
  exit 1
fi

#------------------------------------------------------------------------------
# Auxilary checks. (Assign defaults whenever possible)
#------------------------------------------------------------------------------

if [ -z "$EVALUATION_CONCURRENCY" ]; then
  printf "*** Warning: Evaluation concurrency needs to be defined in the "
  printf "configuration file as well. Should match dakota input file for "
  printf "asynchronous evaluations.\n"
  EVALUATION_CONCURRENCY=1
fi
