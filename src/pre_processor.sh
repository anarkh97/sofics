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

# NOTE: Not strictly necessary, but it's conventional to place 
# structure files in a separate directory for better organization.
struct_dir=$WORKING_DIR/StructModel
struct_inp=$WORKING_DIR/$AEROS_INPUT

gmsh_name=$(basename "$WORKING_DIR/$GMSH_INPUT" .geo)
gmsh_out="${gmsh_name}.msh"

# local log files
gmsh2aeros_log=$struct_dir/log.out

# add evaluation number
sed -i "s/{eval_num}/${DAK_EVAL_NUM}/g" $struct_inp

# make the struct model directory
if [ ! -d "${struct_dir}" ]
then
  mkdir "${struct_dir}"
fi


#------------------------------------------------------------------------------
# Update gmsh input file
#------------------------------------------------------------------------------
substitution_error=0
if [[ ${#cdv[@]} == $num_cdv ]]; then
  for i in "${!cdv[@]}"
  do
    if ! sed -i "s/{cdv_$((i+1))}/${cdv[$i]}/g" "$WORKING_DIR/$GMSH_INPUT"; then
      printf "*** Error: Substitution for continuous design variable cdv_${i} "
      printf "in file ${WORKING_DIR}/${GMSH_INPUT} failed.\n"
      substitution_error=$((substitution_error+1))
    fi
  done
else
  printf "*** Error: Required number of design variables not found.\n"
  exit 1
fi

if [[ ${#csv[@]} == $num_csv ]]; then
  for i in "${!csv[@]}"
  do
    if ! sed -i "s/{csv_$((i+1))}/${csv[$i]}/g" "$WORKING_DIR/$GMSH_INPUT"; then
      printf "*** Error: Substitution for continuous state variable csv_${i} "
      printf "in file ${WORKING_DIR}/${GMSH_INPUT} failed.\n"
      substitution_error=$((substitution_error+1))
    fi
  done
else
  printf "*** Error: Required number of state variables not found.\n"
  exit 1
fi

if [[ $substitution_error -gt 0 ]]; then
  exit 1
fi

#------------------------------------------------------------------------------
# Generate mesh
#------------------------------------------------------------------------------
{
  $GMSH_EXE -3 -format msh -o $struct_dir/$gmsh_out \
    $WORKING_DIR/$GMSH_INPUT > $gmsh2aeros_log
} || {
  # catch error
  printf "*** Error: Failed to generate a mesh.\n"
  exit 1
}

if grep -q "Error" $gmsh2aeros_log; then
  printf "*** Error: GMSH failed to generate the mesh for design "
  printf "${DAK_EVAL_NUM}. Aborting ...\n"
  exit 1
fi 

# convert to aero-s files
$DRIVER_DIR/gmsh2aeros $struct_dir/$gmsh_out $struct_dir/mesh.include >> \
  $gmsh2aeros_log

if grep -q "Error" $gmsh2aeros_log; then
  printf "*** Error: Unable to convert GMSH mesh file to Aero-S inputs for "
  printf "design ${DAK_EVAL_NUM}. Aborting ...\n"
  exit 1
fi 

# get volume elements -- gmsh2aeros only writes volume elements
# TODO: make gmsh2aeros output element number range for general
# representation.
elem_start=1
elem_end=$(
  grep "3D elements" $gmsh2aeros_log |
  cut -d : -f 2 |
  awk '{print $1}'
)

sed -i "s/{elem_start}/${elem_start}/g" $struct_inp
sed -i "s/{elem_end}/${elem_end}/g" $struct_inp

# before launch check and clean-up previous result files
if [ ! -d "$WORKING_DIR/results" ]
then
  mkdir $WORKING_DIR/results
else
  rm -f $WORKING_DIR/results/*
fi

if [ ! -d "$WORKING_DIR/restart" ]
then
  mkdir $WORKING_DIR/restart
else
  rm -f $WORKING_DIR/restart/*
fi
