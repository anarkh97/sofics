#!/bin/bash

STRUCT_DIR=$WORKING_DIR/StructModel
STRUCT_OUT_DIR=$WORKING_DIR/results
POSTPRO=$DRIVER_DIR/postprocessor

# compute mass
$AEROS_EXE $AEROS_INPUT -t

# convert to function values and write to dakota results
# NOTE: `postprocessor` does not know the order in which response
# values must be written. User can call it multiple times, in the
# order that dakota expects and instead append to the file by removing
# the "-w" option in the command line.
$POSTPRO -w -s $STRUCT_OUT_DIR/mass.out -d $DAK_RESULTS \
  > $WORKING_DIR/post_pro_log.out
$POSTPRO -a $STRUCT_OUT_DIR/epstrain -d $DAK_RESULTS \
  >> $WORKING_DIR/post_pro_log.out

if grep -q "Error" $WORKING_DIR/post_pro_log.out; then
  exit 1 
fi
