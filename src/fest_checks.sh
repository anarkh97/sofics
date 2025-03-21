#!/bin/bash

#------------------------------------------------------------------------------
# Only check approximate fluid solver inputs if explicitly requested.
#------------------------------------------------------------------------------

BUILD_DIR=$(dirname "$DRIVER_DIR")
PACKAGE_DIR=$BUILD_DIR/packages

if [ -z "$FEST_SIZE" ]; then
  FEST_SIZE=1
fi

if [ -z "$FEST_EXE" ]; then
  printf "*** Warning: Executable for FEST was "
  printf "not provided. Using default one ...\n"
  if [[ ! -x "$PACKAGE_DIR/fest/fest" ]]; then
    printf "*** Error: Could not find a valid executable for FEST. "
    printf "Aborting ...\n"
    exit 1
  fi
  FEST_EXE="$PACKAGE_DIR/fest/fest"
elif [ -n "$FEST_EXE" ] && [[ ! -x "$FEST_EXE" ]]; then
  printf "*** Error: \"%s\" is not a valid executable. " "${FEST_EXE}"
  printf "Aborting ...\n"
  exit 1
fi

if [[ ! -f "$TEMPLATE_DIR/${FEST_INPUT}.template" ]]; then
  printf "*** Error: Could not find a template file for FEST input "
  printf "file (or template extension is missing).\n"
  exit 1
fi

if [[ ! -f "$TEMPLATE_DIR/${FEST_ERROR_INPUT}.template" ]]; then
  printf "*** Error: Could not find a template file for FEST input "
  printf "file when executed in error mode. (or template extension "
  printf "is missing).\n"
  exit 1
fi
