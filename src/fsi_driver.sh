#!/bin/bash

#------------------------------------------------------------------------------
# Calculate node (list) for the current evaluation.
#------------------------------------------------------------------------------

total_proc=$((
  "$M2C_SIZE"+"$AEROS_SIZE"
))
host_list=$(
  "$DRIVER_DIR/node_list_calculator.sh" "$AEROS_SIZE" "$M2C_SIZE"
)

if [ -z "$host_list" ]; then
  printf "*** Error: Could not calculate the list of host nodes "
  printf "for evaluation %s.\n" "${DAK_EVAL_NUM}"
  exit 1
fi

# -----------------------------------------------------------------------------
# Run analysis in background
# -----------------------------------------------------------------------------
mpiexec --bind-to none \
	-n "$M2C_SIZE" \
  --host "$host_list":"$total_proc" \
	"$M2C_EXE" "$WORKING_DIR/$M2C_INPUT" : \
  -n "$AEROS_SIZE" \
	--host "$host_list":"$total_proc" \
  "$AEROS_EXE" "$WORKING_DIR/$AEROS_INPUT" 2>&1 \
	| tee "$WORKING_DIR/log.out" > /dev/null &

mpi_pid=$!

printf "\033[34mLaunching Evaluation %s on nodes " "${DAK_EVAL_NUM}"
printf "%s with process id %s.\033[0m\n" "${host_list[*]}" "$mpi_pid"

# -----------------------------------------------------------------------------
# Wait for completion
# -----------------------------------------------------------------------------

wait "$mpi_pid"

# -----------------------------------------------------------------------------
# Exit assuming everything was successfull
# -----------------------------------------------------------------------------
exit 0
