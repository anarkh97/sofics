#!/bin/bash

STRUCT_DIR=$WORKING_DIR/StructModel
STRUCT_OUT_DIR=$WORKING_DIR/results
POSTPRO=$DRIVER_DIR/postprocessor

# compute mass
$AEROS_EXE $AEROS_INPUT -t

# convert to function values and write to dakota results
# to append to the results file remove the "-w" option.
$POSTPRO -w -s $STRUCT_OUT_DIR/mass.out -a $STRUCT_OUT_DIR/epstrain \
  -d $DAK_RESULTS > $WORKING_DIR/post_pro_log.out

if grep -q "Error" $WORKING_DIR/post_pro_log.out; then
  exit 1 
fi
