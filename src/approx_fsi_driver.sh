#!/bin/bash

#------------------------------------------------------------------------------
# Calculate node (list) for the current evaluation.
#------------------------------------------------------------------------------

total_proc=$((
  "$FEST_SIZE"+"$AEROS_SIZE"
))
host_list=$(
  "$DRIVER_DIR/node_list_calculator.sh" "$AEROS_SIZE" "$FEST_SIZE"
)
exit_status=$?

# -----------------------------------------------------------------------------
# Run analysis in background
# -----------------------------------------------------------------------------
if [ "$exit_status" -ne 0 ]; then
  printf "*** Error: Could not start the evaluation %s.\n" "${DAK_EVAL_NUM}"
  exit 1
elif [ -n "$host_list" ]; then
  mpiexec --bind-to none \
  	-n "$FEST_SIZE" \
    --host "$host_list":"$total_proc" \
  	"$FEST_EXE" "$WORKING_DIR/$FEST_INPUT" : \
    -n "$AEROS_SIZE" \
  	--host "$host_list":"$total_proc" \
    "$AEROS_EXE" "$WORKING_DIR/$AEROS_INPUT" 2>&1 \
    > >(tee "$WORKING_DIR/log.out" > /dev/null) 2>&1 &
  
  mpi_pid=$!
else
  mpiexec --bind-to none \
  	-n "$FEST_SIZE" \
  	"$FEST_EXE" "$WORKING_DIR/$FEST_INPUT" : \
    -n "$AEROS_SIZE" \
    "$AEROS_EXE" "$WORKING_DIR/$AEROS_INPUT" 2>&1 \
    > >(tee "$WORKING_DIR/log.out" > /dev/null) 2>&1 &
  
  mpi_pid=$!
fi

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
