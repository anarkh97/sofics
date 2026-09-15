#!/bin/bash

struct_out_dir="$WORKING_DIR/results"
postpro="$SOFICS_BIN/aeros2dakota"

# convert to function values and write to dakota results
# NOTE: `postprocessor` does not know the order in which response
# values must be written. User can call it multiple times, in the
# order that dakota expects and instead append to the file by removing
# the "-w" option in the command line.
"$postpro" -w -v "$struct_out_dir/velo.0" -d "$DAK_RESULTS" \
  > "$WORKING_DIR/post_pro_log.out"
"$postpro" -v "$struct_out_dir/disp.0" -d "$DAK_RESULTS" \
  >> "$WORKING_DIR/post_pro_log.out"

if grep -q "Error" "$WORKING_DIR/post_pro_log.out"; then
  exit 1 
fi

# manualy scaling objective and constraint values
readarray -t values < <(awk '{print $1}' "$DAK_RESULTS")

objective=$(awk "BEGIN {print ${values[0]} / 30000.0}")
constraint=$(awk "BEGIN {print ${values[1]} / 1.0}")

printf "    %.6e\n    %.6e\n" "$objective" \
  "$constraint" > "$DAK_RESULTS"

exit 0
