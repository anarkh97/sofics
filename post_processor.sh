#!/bin/bash

STRUCT_DIR=$WORKING_DIR/StructModel
STRUCT_OUT_DIR=$WORKING_DIR/results
POSTPRO=/projects/wang_aoe_lab/AdityaNarkhede/DakotaOptimization/CodeRepository

# compute mass
$AEROS_EXE $AEROS_INPUT -t

# compute max plastic strain and output results
$POSTPRO/PostProcessor/post_processor \
  $STRUCT_DIR/mesh.include.surf6 \
  $STRUCT_OUT_DIR/mass.out $STRUCT_OUT_DIR/epstrain \
  $DAK_RESULTS
