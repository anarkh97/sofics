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

# NOTE: Not strictly necessary, but it's conventional to place 
# structure files in a separate directory for better organization.
struct_dir=$WORKING_DIR/StructModel
struct_inp=$WORKING_DIR/$AEROS_INPUT

# make the struct model directory
if [ ! -d "${struct_dir}" ]
then
  mkdir "${struct_dir}"
fi

#------------------------------------------------------------------------------
# Copy mesh
#------------------------------------------------------------------------------

cp "$TEMPLATE_DIR/mesh.include" "$struct_dir/mesh.include"
cp "$TEMPLATE_DIR/mesh.include.surf1" "$struct_dir/mesh.include.surf1"

#------------------------------------------------------------------------------
# Update Aero-S input file.
#------------------------------------------------------------------------------

# add evaluation number
sed -i "s/{eval_num}/${DAK_EVAL_NUM}/g" "$struct_inp"

{
  "$LAUNCH_DIR"/attribute_generator "$struct_dir/mesh.include" \
    "$struct_dir/matlaw.include" "${cdv[@]}" > "$struct_dir/attribute_gen.out"
} || {
  printf "*** Error: Failed while generating material attributes.\n"
  exit 1
}

#------------------------------------------------------------------------------
# Before launch check and clean-up previous result files
#------------------------------------------------------------------------------
if [ ! -d "$WORKING_DIR/results" ]
then
  mkdir "$WORKING_DIR/results"
else
  rm -f "$WORKING_DIR/results/*"
fi

if [ ! -d "$WORKING_DIR/restart" ]
then
  mkdir "$WORKING_DIR/restart"
else
  rm -f "$WORKING_DIR/restart/*"
fi

#------------------------------------------------------------------------------
# Exit
#------------------------------------------------------------------------------
exit 0
