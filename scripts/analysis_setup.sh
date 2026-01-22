#!/bin/bash

#------------------------------------------------------------------------------
# Perform checks before proceeding (in CURRENT shell)
#------------------------------------------------------------------------------
if ! source "${SOFICS_BIN}/fsi_checks.sh"; then
  # NOTE: Here we do not write a results file and force dakota to fail.
  # This is a brute force approach for cases when user specifies a failure
  # capture method other than `abort` in dakota.

  exit 1
fi

#------------------------------------------------------------------------------
# Internal variables (local to current shell)
#------------------------------------------------------------------------------
neighbors=""
target_id=""

#------------------------------------------------------------------------------
# Read Dakota Parameter file and forward variables to a pre-processor
#------------------------------------------------------------------------------
{ 
  num_var=$(awk 'NR==1 {print $1}' "$WORKING_DIR/$DAK_PARAMS")
} || {
  printf "*** Error: Something went wrong while extracting number "
  printf "of variables from \"%s\".\n" "$WORKING_DIR/$DAK_PARAMS"
  exit 1
}

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
cdv_list=""
for i in $(seq 1 "$num_cdv")
do
  val=$(awk -v line="$i" 'NR==line+1 {print $1}' "$WORKING_DIR/$DAK_PARAMS")
  desc=$(awk -v line="$i" 'NR==line+1 {print $2}' "$WORKING_DIR/$DAK_PARAMS")

  # This is to ensure that only continuous design variables type were
  # used in the Dakota input file. Will be extended in future and this condition
  # will be relaxed to accommodate discrete variables of type INTEGER and STRING. 
  if [[ "$desc" != "cdv_$i" ]]; then
    printf "*** Error: This utility expects continuous design "
    printf "variables with dakota's default descriptors. Found a variable "
    printf "with descriptor %s. Aborting ...\n" "$desc"
    exit 1
  fi

  if [ -z "$cdv_list" ]; then
    cdv_list="$val"
  else
    cdv_list="$cdv_list:$val"
  fi
done

# These are state variables that the user can use to define psuedo/problem 
# specific values, such as mesh size, directly from Dakota's input. We expect 
# these to be of type REAL, i.e. continuous state variables, with an expections
# for a special (internal) state variable with KEYWORDS: "SWITCH", "NEIGHBOR".
num_sv=$((num_var-num_cdv))
num_csv=0
csv_list=""
for i in $(seq 1 $num_sv)
do
  offset=$((num_cdv+i))
  val=$(awk -v line=$offset 'NR==line+1 {print $1}' "$WORKING_DIR/$DAK_PARAMS")
  desc=$(awk -v line=$offset 'NR==line+1 {print $2}' "$WORKING_DIR/$DAK_PARAMS")

  # This is to ensure that only continuous state variables of REAL type were
  # used in the Dakota input file. We internally handle INTEGER and STRING variables
  # used for ADOPT.
  if [[ $desc != "csv_$i" ]]; then

    # SWITCH solver flag for ADOPT.
    if [[ $desc == "SWITCH" ]]; then
      # check for our special keywords
      if [[ $val == "TRUE" || $val == "APPROX" || $val == "ERROR" ]]; then
        SOLVER_TYPE="$val"
      else
        printf "*** Error: The SWTICH state variable is a reserved keyword " 
        printf "with \"TRUE\", \"APPROX\", and \"ERROR\" options. Instead "
        printf "recieved %s. Aborting ...\n" "$val"
        exit 1
      fi
      continue
    fi

    # Integer state used to pass the nearest neighbors (regex match).
    if [[ $desc =~ ^NEIGHBOR_[0-9]+$ ]]; then
     
      # (AN): The following in no longer valid. ADOPT uses discrete_state_range
      # ADOPT uses continuous state variables to handle evaluation IDs.
      # These variables are of type "REAL" and are written in scientific
      # format. Here were correct them to integers.
      val=$(awk -v value="$val" 'BEGIN {print int(value)}')

      # assign if empty
      if [ -z "$neighbors" ]; then
        neighbors="$val"
      else
        neighbors="$neighbors:$val" # append
      fi

      continue
    fi

    # Target evaluation id --- used when performing error simulations.
    if [[ "$desc" == "TARGET" ]]; then

      val=$(awk -v value="$val" 'BEGIN {print int(value)}')
      # assign if empty
      if [ -z "$target_id" ]; then
        target_id="$val"
      else
        printf "*** Error: Found multiple targets in the parameter file.\n"
        exit 1
      fi

      continue
    fi

    # Other types are "errors" for SOFICS
    printf "*** Error: This utility expects continuous state "
    printf "variables with dakota's default descriptors. Found a variable "
    printf "with descriptor %s. Aborting ...\n" "$desc"
    exit 1

  fi

  if [ -z "$csv_list" ]; then
    csv_list="$val"
  else
    csv_list="$csv_list:$val"
  fi
  num_csv=$((num_csv+1))

done

#------------------------------------------------------------------------------
# Perform FEST related checks and generate "META" file
#------------------------------------------------------------------------------
if [[ "$SOLVER_TYPE" == "APPROX" || "$SOLVER_TYPE" == "ERROR" ]]; then

  # Run checks in CURRENT shell.
  if ! source "${SOFICS_BIN}/approx_fsi_checks.sh"; then
    # NOTE: Here we do not write a results file and force dakota to fail.
    # This is a brute force approach for cases when user specifies a failure
    # capture method other than `abort` in dakota.
  
    exit 1
  fi

  # This will be created in the WORKING_DIR
  if ! bash "$META_GENERATOR_FILE" "$target_id" "$cdv_list" "$neighbors"
  then
    printf "*** Error: Failed to generate metafile for "
    printf "ADOPT (design %s).\n" "${DAK_EVAL_NUM}"

    printf "FAIL\n" > "$DAK_RESULTS"
    exit 0
  fi

fi

#------------------------------------------------------------------------------
# Copy all templates to the working directory
#------------------------------------------------------------------------------
if [[ "$SOLVER_TYPE" != "ERROR" ]]; then
  cp "$TEMPLATE_DIR/${AEROS_INPUT}.template" "$WORKING_DIR/$AEROS_INPUT"
fi

# copy the fluid files based on solver type
case "$SOLVER_TYPE" in

  TRUE)
    cp "$TEMPLATE_DIR/${M2C_INPUT}.template" "$WORKING_DIR/$M2C_INPUT"
    # copy any auxilary inputs that were provided for m2c
    if [ -n "$M2C_AUX" ]; then
      # get names of all auxilary inputs
      IFS=: read -r -a fluid_aux_inps <<< "$M2C_AUX"
    
      # remove empty fields (or spaces)
      fluid_aux_inps=("${fluid_aux_inps[@]// /}")
    
      for i in "${!fluid_aux_inps[@]}"; do
        cp "$TEMPLATE_DIR/${fluid_aux_inps[$i]}" \
          "$WORKING_DIR/${fluid_aux_inps[$i]}"
      done
    fi
    ;;
  APPROX)
    cp "$TEMPLATE_DIR/${FEST_INPUT}.template" "$WORKING_DIR/$FEST_INPUT"
    ;;
  ERROR)
    cp "$TEMPLATE_DIR/${FEST_ERROR_INPUT}.template" \
      "$WORKING_DIR/$FEST_ERROR_INPUT"
    ;;

esac

#------------------------------------------------------------------------------
# Execute structural pre-processor in a sub-shell
#------------------------------------------------------------------------------

# Error run does not require structural mesh.
if [[ "$SOLVER_TYPE" != "ERROR" ]]; then
  if ! bash "$PREPROCESS_FILE" "$cdv_list" "$csv_list" ; then
    printf "*** Error: Failed at pre-processing stage for design "
    printf "%s.\n" "${DAK_EVAL_NUM}"
  
    # Here we let dakota capture the failure and proceed 
    # based on user specification.
  
    printf "FAIL\n" > "$DAK_RESULTS"
    exit 0
  fi
fi
