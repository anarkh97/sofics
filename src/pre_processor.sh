#!/bin/bash

# All global variables are depicted in capital text. For a list
# of global variables check the `driver.sh` file.

# NOTE: Not strictly necessary, but it's conventional to place 
# structure files in a separate directory for better organization.
STRUCT_DIR=$WORKING_DIR/StructModel
STRUCT_INP=$WORKING_DIR/fem.in

# local log files
gmsh2aeros_log=$STRUCT_DIR/log.out

# add evaluation number
sed -i "s/{eval_num}/${DAK_EVAL_NUM}/g" $WORKING_DIR/fem.in

# make the struct model directory
if [ ! -d "${STRUCT_DIR}" ]
then
  mkdir "${STRUCT_DIR}"
fi

# parse the dakota params file
num_dv=$(awk 'NR==1 {print $1}' "$WORKING_DIR/$DAK_PARAMS")
cdv=()
for i in $(seq 1 $num_dv)
do
  val=$(awk -v line=$i 'NR==line+1 {print $1}' "$WORKING_DIR/$DAK_PARAMS")
  desc=$(awk -v line=$i 'NR==line+1 {print $2}' "$WORKING_DIR/$DAK_PARAMS")

  if [[ $desc != "cdv_$i" ]]; then
    echo "*** Error: This utility expects the user to use Dakota's default \
      continuous design variable descriptors. Remove any special descriptors \
      and try again. Aborting ..."
    exit 1
  fi

  cdv+=($val)
done

# update gmsh input file
if [[ "${#cdv[@]}" == $num_dv ]]
then
  for i in $(seq 1 $num_dv)
  do
    sed -i "s/{cdv_$i}/${cdv[$i]}/g" $WORKING_DIR/struct.geo
  done
else
  echo "***Error: Required number of design variables not found."
  exit 1
fi

# generate mesh
gmsh -3 -format msh -o $STRUCT_DIR/struct.msh $WORKING_DIR/struct.geo > \
  $gmsh2aeros_log

if grep -q "Error" $gmsh2aeros_log; then
  echo "*** Error: GMSH failed to generate the mesh for design \
    ${DAK_EVAL_NUM}. Aborting ..."
  exit 1
fi 

# convert to aero-s files
gmsh2aeros $STRUCT_DIR/struct.msh $STRUCT_DIR/mesh.include >> \
  $gmsh2aeros_log

if grep -q "Error" $gmsh2aeros_log; then
  echo "*** Error: Unable to convert GMSH mesh file to Aero-S inputs for \
    design ${DAK_EVAL_NUM}. Aborting ..."
  exit 1
fi 

# get volume elements -- gmsh2aeros only writes volume elements
elem_start=$(
  grep "3D elements" $gmsh2aeros_log |
  cut -d : -f 1 |
  awk 'print $1}'
)
elem_end=$(
  grep "3D elements" $gmsh2aeros_log |
  cut -d : -f 2 |
  awk '{print $1}'
)

sed -i "s/{elem_start}/${elem_start}/g" $WORKING_DIR/fem.in
sed -i "s/{elem_end}/${elem_end}/g" $WORKING_DIR/fem.in

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
