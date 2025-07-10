#!/bin/bash

#------------------------------------------------------------------------------
# Calculate node (list) for the current evaluation.
#------------------------------------------------------------------------------

host_list=""
total_proc=1
if ! source "${DRIVER_DIR}/node_list_calculator.sh" "$host_list" "0" "1"
then
  printf "*** Error: Could not calculate the list of host nodes "
  printf "for (error) evaluation %s.\n" "${DAK_EVAL_NUM}"
fi

### run analysis in background
mpiexec --bind-to none \
	-n "$FEST_SIZE" \
  --host "$host_list":"$total_proc" \
	"$FEST_EXE" "$WORKING_DIR/$FEST_ERROR_INPUT" \
	| tee "$WORKING_DIR/log.out" > /dev/null &

mpi_pid=$!

printf "\033[34mLaunching (error) Evaluation %s on nodes " "${DAK_EVAL_NUM}"
printf "%s with process id %s.\033[0m\n" "${host_list[*]}" "$mpi_pid"

# -----------------------------------------------------------------------------
# Wait for completion
# -----------------------------------------------------------------------------

wait "$mpi_pid"

# -----------------------------------------------------------------------------
# Exit assuming everything was successfull
# -----------------------------------------------------------------------------
exit 0
