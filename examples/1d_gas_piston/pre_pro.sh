#!/bin/bash

cont_design_vars=$1
cont_state_vars=$2

#------------------------------------------------------------------------------
# Convert ":" separated lists to an array
#------------------------------------------------------------------------------
if [ -z "$cont_design_vars" ]
then
  printf "*** Error: Empty list of continuous design variables.\n"
  exit 1
else
  IFS=: read -r -a cdv <<< "$cont_design_vars"  

  # remove empty fields (or spaces)
  cdv=("${cdv[@]// /}")
fi

if [ -n "$cont_state_vars" ]
then
  IFS=: read -r -a csv <<< "$cont_state_vars"

  # remove empty fields (or spaces)
  csv=("${csv[@]// /}")
fi

#------------------------------------------------------------------------------
# Sub-directory for structural mesh.
#------------------------------------------------------------------------------

equiv_a=$(
  awk "BEGIN { printf \"%16.8e\", 1/(1/${cdv[0]} + 1/${cdv[1]}) }"
)
equiv_k=$(
  awk "BEGIN { printf \"%.3e\", (1e5)*(${equiv_a}) }"
)
equiv_m=$(
  awk "BEGIN { printf \"%.3e\", (1e-3)*(${cdv[0]}+${cdv[1]}) }"
)

sed -i "s/{equiv_k}/${equiv_k}/g" "$AEROS_INPUT"
sed -i "s/{equiv_m}/${equiv_m}/g" "$AEROS_INPUT"

#------------------------------------------------------------------------------
# Before launch check and clean-up previous result files
#------------------------------------------------------------------------------
if [ ! -d "$WORKING_DIR/results" ]
then
  mkdir "$WORKING_DIR/results"
else
  rm -f "$WORKING_DIR/results/*"
fi

#------------------------------------------------------------------------------
# Exit
#------------------------------------------------------------------------------
exit 0
