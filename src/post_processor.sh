#!/bin/bash

STRUCT_DIR=$WORKING_DIR/StructModel
STRUCT_OUT_DIR=$WORKING_DIR/results
POSTPRO=$DRIVER_DIR/postprocessor

# compute mass
$AEROS_EXE $AEROS_INPUT -t

# convert to function values and write to dakota results
$POSTPRO $STRUCT_DIR/mesh.include.surf6 \
  $STRUCT_OUT_DIR/mass.out $STRUCT_OUT_DIR/epstrain \
  $DAK_RESULTS > $WORKING_DIR/post_pro_log.out

if grep -q "Error" $WORKING_DIR/post_pro_log.out; then
  exit 1 
fi
