#!/bin/bash

#------------------------------------------------------------------------------
# Only check approximate fluid solver inputs if explicitly requested.
#------------------------------------------------------------------------------

if [ -z "$FEST_SIZE" ]; then
  FEST_SIZE=1
fi

if [ -z "$FEST_EXE" ]; then
  printf "*** Warning: Executable for FEST was "
  printf "not provided. Using default one ...\n"
  if [[ ! -x "$DEFAULT_FEST_EXE" ]]; then
    printf "*** Error: Could not find a valid executable for FEST. "
    printf "Aborting ...\n"
    exit 1
  fi
  FEST_EXE="$DEFAULT_FEST_EXE"
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

if [[ -z "$META_FILE" ]]; then
  printf "*** Warning: Name for file containing nearest neighbor "
  printf "information was not provided. Using default "
  printf "\"meta_data.txt\".\n"

  META_FILE="meta_data.txt"
fi

if [[ -z "$META_SURFACE_FILE" ]]; then
  printf "*** Warning: Name for file containing embedded surface "
  printf "coordinates was not provided. Using default "
  printf "\"surface.top\".\n"

  META_SURFACE_FILE="surface.top"
fi

if [[ -z "$META_SOLUTION_FILE" ]]; then
  printf "*** Warning: Name for file containing pressure time "
  printf "history for the embedded surface was not provided. "
  printf "Using default \"surface_load_2.txt\".\n"

  META_SOLUTION_FILE="surface_load_2.txt"
fi
