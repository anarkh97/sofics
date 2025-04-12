#!/bin/bash

#------------------------------------------------------------------------------
# Inputs
#------------------------------------------------------------------------------
targ_eval_id=$1
targ_vars_list=$2
neighbors_list=$3

#------------------------------------------------------------------------------
# Convert to arrays
#------------------------------------------------------------------------------
IFS=: read -r -a targ_vars <<< "$targ_vars_list"
IFS=: read -r -a neighbors <<< "$neighbors_list"

# remove empty fields (or spaces)
targ_vars=("${targ_vars[@]// /}")
neighbors=("${neighbors[@]// /}")

# number of continuous optimization variables
num_cont_vars="${#targ_vars[@]}"

#------------------------------------------------------------------------------
# Generate meta file and headers
#------------------------------------------------------------------------------
if [ -f "$WORKING_DIR/$META_FILE" ]; then
  # clear any exiting meta files.
  rm -rf "$WORKING_DIR/$META_FILE"
fi

# write headers
printf "## Meta file containing paths to data stored for nearest neighbors\n" \
  > "$WORKING_DIR/$META_FILE"
printf "## Type  Directory  Mesh  Solution  Parameters\n" \
  >> "$WORKING_DIR/$META_FILE"

#------------------------------------------------------------------------------
# Write neighbors to metafile.
#------------------------------------------------------------------------------
directory_pattern="${WORKING_DIR##*/}" # with the eval tag
directory_pattern="${directory_pattern%[0-9]*}" # without the eval tag

param_pattern=$(basename "$DAK_PARAMS") # with the eval tag
param_pattern="${param_pattern%[0-9]*}" # without the eval tag

relative_path=$(
  realpath -s --relative-to="$WORKING_DIR" "$LAUNCH_DIR"
)

for n in "${!neighbors[@]}"; do

  # first get parameters -- assuming eval tags.
  neighbor_eval_id="${neighbors[$n]}"
  neighbor_direct="$relative_path/${directory_pattern}${neighbor_eval_id}"
  neighbor_params="$neighbor_direct/${param_pattern}${neighbor_eval_id}"

  # read parameters
  num_var=$(awk 'NR==1 {print $1}' "$neighbor_params")
  if [[ $num_var -lt $num_cont_vars ]]; then
    printf "*** Error: Enough continuous design variables were not provided "
    printf "for evaluation ID %s. Aborting ...\n" "$neighbor_eval_id"
    exit 1
  fi
  
  # we only parse the design variables and ignore the rest.
  # NOTE: dakota always writes the continuous design variables first.
  # TODO: Extend to include other kinds of design variables as well.
  # TODO: Could remove these checks as these have already been done
  #       in analysis setup.
  neighbor_cdv=()
  for i in $(seq 1 "$num_cont_vars")
  do
    val=$(awk -v line="$i" 'NR==line+1 {print $1}' "$neighbor_params")
    desc=$(awk -v line="$i" 'NR==line+1 {print $2}' "$neighbor_params")
  
    # This is to ensure that only continuous design variables type were
    # used in the Dakota input file. Will be extended in future and this condition
    # will be relaxed to accommodate discrete variables of type INTEGER and STRING. 
    if [[ "$desc" != "cdv_$i" ]]; then
      printf "*** Error: This utility expects continuous design "
      printf "variables with dakota's default descriptors. Found a variable "
      printf "with descriptor %s. Aborting ...\n" "$desc"
      exit 1
    fi
 
    neighbor_cdv+=("$val")

  done

  # write this neighbor to metafile
  {
    printf "NEIGHBOR  "
    printf "\"%s/\"  " "$neighbor_direct"
    printf "\"%s\"  " "$META_SURFACE_FILE"
    printf "\"%s\"  " "$META_SOLUTION_FILE"
    # Pipe the output of printf to sed to remove trailing whitespaces
    printf "%s  " "${neighbor_cdv[@]}" | sed 's/[[:space:]]*$//'
    printf "\n"
  } >> "$META_FILE"


done

#------------------------------------------------------------------------------
# Write target to metafile.
#------------------------------------------------------------------------------
{
  printf "TARGET  "
  printf "\"%s/\"  " "$relative_path/${directory_pattern}${targ_eval_id}"
  printf "\"%s\"  " "$META_SURFACE_FILE"
  printf "\"%s\"  " "$META_SOLUTION_FILE"
  printf "%s  " "${targ_vars[@]}" | sed 's/[[:space:]]*$//'
  printf "\n"
} >> $META_FILE

#------------------------------------------------------------------------------
# Exit successfully
#------------------------------------------------------------------------------
exit 0
