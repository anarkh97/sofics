#!/bin/bash

# All global variables are depicted in capital text. For a list
# of global variables check the `driver.sh` file.

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

# parse the dakota params file
num_var=$(awk 'NR==1 {print $1}' "$WORKING_DIR/$DAK_PARAMS")
# Dakota specifies continous design variables as derivative variables
# even if gradient-free optimization is employed. There is no other 
# way to explicitly distinguish design variables from state variables.
num_cdv=$(
  grep "derivative_variables" "$WORKING_DIR/$DAK_PARAMS" |
  awk '{print $1}'
)

if [[ $num_var -lt $num_cdv ]]; then
  printf "*** Error: Enough continuous design variables were not provided. "
  printf "Aborting ...\n"
  exit 1
fi

# we only parse the design variables and ignore the rest.
# NOTE: dakota always writes the continuous design variables first.
# TODO: Extend to include other kinds of design variables as well.
cdv=()
for i in $(seq 1 $num_cdv)
do
  val=$(awk -v line=$i 'NR==line+1 {print $1}' "$WORKING_DIR/$DAK_PARAMS")
  desc=$(awk -v line=$i 'NR==line+1 {print $2}' "$WORKING_DIR/$DAK_PARAMS")

  # This is to ensure that only continuous design variables type were
  # used in the Dakota input file. Will be extended in future and this condition
  # will be relaxed to accommodate discrete variables of type INTEGER and STRING. 
  if [[ $desc != "cdv_$i" ]]; then
    printf "*** Error: This utility expects continuous design "
    printf "variables with dakota's default descriptors. Found a variable "
    printf "with descriptor %s. Aborting ...\n" "$desc"
    exit 1
  fi

  cdv+=($val)
done

# These are state variables that the user can use to define psuedo/problem 
# specific values, such as mesh size, directly from Dakota's input. We still
# expect these to be of type REAL, i.e. continuous state variables.
num_csv=$((num_var-num_cdv))
csv=()
for i in $(seq 1 $num_csv)
do
  offset=$((num_cdv+i))
  val=$(awk -v line=$offset 'NR==line+1 {print $1}' "$WORKING_DIR/$DAK_PARAMS")
  desc=$(awk -v line=$offset 'NR==line+1 {print $2}' "$WORKING_DIR/$DAK_PARAMS")

  # This is to ensure that only continuous state variables of REAL type were
  # used in the Dakota input file. Will be extended in future and this condition
  # will be relaxed to accommodate discrete variables of type INTEGER and STRING. 
  if [[ $desc != "csv_$i" ]]; then
    printf "*** Error: This utility expects continuous state "
    printf "variables with dakota's default descriptors. Found a variable "
    printf "with descriptor %s. Aborting ...\n" "$desc"
    exit 1
  fi

  csv+=($val)
done

# update gmsh input file
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

# generate mesh
$GMSH_EXE -3 -format msh -o $struct_dir/$gmsh_out $WORKING_DIR/$GMSH_INPUT > \
  $gmsh2aeros_log

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
