#!/bin/bash

# All global variables are depicted in capital text. For a list
# of global variables check the `driver.sh` file.

# NOTE: Not strictly necessary, but it's conventional to place 
# structure files in a separate directory for better organization.
STRUCT_DIR=$WORKING_DIR/StructModel
STRUCT_INP=$WORKING_DIR/fem.in

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
  cdv+=($val)
done

# update gmsh input file
if [ "${#cdv[@]}" -ge 2 ]
then
  sed -i "s/{lx}/${cdv[0]}/g" $WORKING_DIR/struct.geo
  sed -i "s/{ly}/${cdv[1]}/g" $WORKING_DIR/struct.geo
else
  echo "***Error: Required number of design variables not found."
  exit 1
fi

# generate mesh
gmsh -3 -format msh -o $STRUCT_DIR/struct.msh $WORKING_DIR/struct.geo 

# convert to aero-s files
gmsh2aeros $STRUCT_DIR/struct.msh $STRUCT_DIR/mesh.include > \
  $STRUCT_DIR/gmsh_log.out

# get volume elements -- gmsh2aeros only writes volume elements
elem_start=1
elem_end=$(
  grep "3D elements" $STRUCT_DIR/gmsh_log.out |
  cut -d : -f 2 |
  awk '{print $1}'
)

sed -i "s/{elem_start}/${elem_start}/g" $WORKING_DIR/fem.in
sed -i "s/{elem_end}/${elem_end}/g" $WORKING_DIR/fem.in

# before launch check and clean-up previous result files
if [ ! -d "${WORKING_DIR}/results" ]
then
  mkdir $WORKING_DIR/results
else
  rm -f $WORKING_DIR/results/*
fi

if [ ! -d "${WORKING_DIR}/restart" ]
then
  mkdir $WORKING_DIR/restart
else
  rm -f $WORKING_DIR/restart/*
fi
