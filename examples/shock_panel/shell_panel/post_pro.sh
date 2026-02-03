#!/bin/bash

struct_out_dir="$WORKING_DIR/results"
postpro="$DRIVER_DIR/postprocessor"

# compute mass
"$AEROS_EXE" "$AEROS_INPUT" -t

"$postpro" -w -s "$struct_out_dir/mass.out" -d "$DAK_RESULTS" \
  > "$WORKING_DIR/post_pro_log.out"

if grep -q "Error" "$WORKING_DIR/post_pro_log.out"; then
  exit 1 
fi

# manualy scaling objective and constraint values
readarray -t values < <(awk '{print $1}' "$DAK_RESULTS")

printf "    %.6e\n" "$response" > "$DAK_RESULTS"

exit 0
