#!/bin/bash

# -----------------------------------------------------------------------------
# SOFICS Build Script
# -----------------------------------------------------------------------------

# Default values
CURRENT_DIR=$(pwd)
BUILD_TYPE=${1:-Release}
INSTALL_PREFIX=${2:-$CURRENT_DIR/sofics}
BUILD_AEROS=${3:-ON}
BUILD_M2C=${4:-ON}
BUILD_TOOLS=${5:-ON}
UPDATE_SUBMODULES=${6:-OFF}
NUM_CPU=${3:-$(nproc)}

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NO_COLOR='\033[0m' # No Color

echo -e "${GREEN}===================${NO_COLOR}"
echo -e "${GREEN}SOFICS Build Script${NO_COLOR}"
echo -e "${GREEN}===================${NO_COLOR}"
echo "Build type:      ${BUILD_TYPE}"
echo "Install prefix:  ${INSTALL_PREFIX}"
echo ""

# Configure
echo -e "${YELLOW} Configuring with CMake ...${NO_COLOR}"
mkdir -p build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
    -DBUILD_AEROS=$BUILD_AEROS \
    -DBUILD_M2C=$BUILD_M2C \
    -DBUILD_TOOLS=$BUILD_TOOLS \
    -DUPDATE_SUBMODULES=$UPDATE_SUBMODULES
echo -e "${GREEN} ... Done.${NO_COLOR}"
echo ""

cd ..

# Install
echo -e "${YELLOW} Installing...${NO_COLOR}"
cmake --install build
echo -e "${GREEN} ... Done.${NO_COLOR}"
echo ""

# Summary
echo -e "${GREEN}========================================${NO_COLOR}"
echo -e "${GREEN}Build Complete!${NO_COLOR}"
echo -e "${GREEN}========================================${NO_COLOR}"
echo ""
echo "Installation directory: ${INSTALL_PREFIX}"
echo ""
echo "Available executables:"
if [ "$BUILD_AEROS" -eq "ON" ]; then
  echo "  - aeros"
fi
if [ "$BUILD_M2C" -eq "ON" ]; then
  echo "  - m2c"
fi
if [ "$BUILD_TOOLS" -eq "ON" ]; then
  echo "  - gmsh2aeros"
  echo "  - aeros2dakota"
fi
echo -e "${GREEN}========================================${NO_COLOR}"
