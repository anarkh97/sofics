#!/bin/bash

struct_out_dir="$WORKING_DIR/results"
postpro="$SOFICS_BIN/postprocessor"

# compute mass
"$AEROS_EXE" "$AEROS_INPUT" -t

# convert to function values and write to dakota results
# NOTE: `postprocessor` does not know the order in which response
# values must be written. User can call it multiple times, in the
# order that dakota expects and instead append to the file by removing
# the "-w" option in the command line.
"$postpro" -w -v "$struct_out_dir/disp.3" -d "$DAK_RESULTS" \
  > "$WORKING_DIR/post_pro_log.out"
"$postpro" -p "$struct_out_dir/vmstress.1" -d "$DAK_RESULTS" \
  >> "$WORKING_DIR/post_pro_log.out"
"$postpro" -s "$struct_out_dir/mass.out" -d "$DAK_RESULTS" \
  >> "$WORKING_DIR/post_pro_log.out"

if grep -q "Error" "$WORKING_DIR/post_pro_log.out"; then
  exit 1 
fi

# manualy scaling objective and constraint values
readarray -t values < <(awk '{print $1}' "$DAK_RESULTS")

objective=$(awk "BEGIN {print ${values[0]} / 6.0}")
constraint1=$(awk "BEGIN {print ${values[1]} / 8.0e8}")
constraint2=$(awk "BEGIN {print ${values[2]} / 0.76}")

printf "    %.6e\n    %.6e\n    %.6e\n" "$objective" \
  "$constraint1" "$constraint2" > "$DAK_RESULTS"

exit 0
