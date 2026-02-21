#!/bin/bash

# -----------------------------------------------------------------------------
# Build defaults and options
# -----------------------------------------------------------------------------

usage() {
  echo "Usage: ./build.sh [options]"
  echo "  --build-type         Build type (default is Release)"
  echo "  --prefix             Install prefix (default is current-dir/build/sofics)"
  echo "  --no-aeros           Disable BUILD_AEROS"
  echo "  --no-m2c             Disable BUILD_M2C"
  echo "  --no-tools           Disable BUILD_TOOLS"
  echo "  --update-submodules  Enable UPDATE_SUBMODULES"
  echo "  --jobs               Number of CPUs (default: 4)"
  echo "  --help               Show this help message"
  exit 1
}

# Default values
CURRENT_DIR=$(pwd)
BUILD_TYPE="Release"
BUILD_DIR=""
INSTALL_PREFIX="$CURRENT_DIR/build/sofics"
BUILD_AEROS=ON
BUILD_M2C=ON
BUILD_TOOLS=ON
UPDATE_SUBMODULES=OFF
NUM_CPU=4

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NO_COLOR='\033[0m' # No Color

# parse command-line arguments using "getopt"
SHORT_OPTIONS=""
LONG_OPTIONS="help,build-type:,build-dir:,prefix:,no-aeros,no-m2c,no-tools,update-submodules,jobs:"

PARSER="$(getopt -o "$SHORT_OPTIONS" -l "$LONG_OPTIONS" -- "$@")" || {
  # getopt failed
  printf "*** Error: Could not understand the option \"%s\".\n" "$?"
  exit 1
}

eval set -- "$PARSER"

# Get user specified options.
while true; do

  case "$1" in

    --help)
      usage
      ;;
    --build-type)
      BUILD_TYPE="$2"
      shift 2
      ;;
    --build-dir)
      BUILD_DIR="$2"
      shift 2
      ;;
    --prefix)
      INSTALL_PREFIX="$2"
      shift 2
      ;;
    --no-aeros)
      BUILD_AEROS=OFF
      ;;
    --no-m2c)
      BUILD_M2C=OFF
      ;;
    --no-tools)
      BUILD_TOOLS=OFF
      ;;
    --update-submodules)
      UPDATE_SUBMODULES=ON
      ;;
    --jobs)
      NUM_CPU="$2"
      shift 2
      ;;
    --)
      shift
      break
      ;;
    *)
      usage 
      ;;

  esac

done

# -----------------------------------------------------------------------------
# Start the build process 
# -----------------------------------------------------------------------------

echo -e "${GREEN}===================${NO_COLOR}"
echo -e "${GREEN}SOFICS Build Script${NO_COLOR}"
echo -e "${GREEN}===================${NO_COLOR}"
echo "Build type:      ${BUILD_TYPE}"
echo "Install prefix:  ${INSTALL_PREFIX}"
echo ""

if [[ -z "$BUILD_DIR" ]]; then
  BUILD_DIR=$(dirname "$INSTALL_PREFIX")
fi

# -----------------------------------------------------------------------------
# Configure (skip if previous cache found)
# -----------------------------------------------------------------------------

if [[ -e "$BUILD_DIR/CMakeCache.txt" ]]; then
  echo -e "${GREEN}Found previous configuration, skipping CMake...${NO_COLOR}"
else
  echo -e "${YELLOW}Configuring with CMake...${NO_COLOR}"
  mkdir -p "$BUILD_DIR" && cd "$BUILD_DIR" || exit 1
  cmake .. \
    -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DBUILD_AEROS="${BUILD_AEROS}" \
    -DBUILD_M2C="${BUILD_M2C}" \
    -DBUILD_TOOLS="${BUILD_TOOLS}" \
    -DUPDATE_SUBMODULES="${UPDATE_SUBMODULES}" \
    || exit 1
  echo -e "${GREEN}Done.${NO_COLOR}\n"
fi

cd "$BUILD_DIR" || exit 1

# -----------------------------------------------------------------------------
# Build targets (parallel)
# -----------------------------------------------------------------------------

echo -e "${YELLOW}Building with $NUM_CPU CPUs...${NO_COLOR}"
make -j"$NUM_CPU" || exit 1
echo -e "${GREEN}Done.${NO_COLOR}\n"

# -----------------------------------------------------------------------------
# Install targets into the install prefix directory
# -----------------------------------------------------------------------------

echo -e "${YELLOW}Installing...${NO_COLOR}"
make install || exit 1
echo -e "${GREEN}Done.${NO_COLOR}\n"

# -----------------------------------------------------------------------------
# Print final build summary
# -----------------------------------------------------------------------------

echo -e "${GREEN}========================================${NO_COLOR}"
echo -e "${GREEN}Build Complete!${NO_COLOR}"
echo -e "${GREEN}========================================${NO_COLOR}"
echo ""
echo "Executables located in: ${INSTALL_PREFIX}/bin"
echo ""
echo "Available executables:"
if [[ "$BUILD_AEROS" -eq "ON" ]]; then
  echo "  - aeros"
fi
if [[ "$BUILD_M2C" -eq "ON" ]]; then
  echo "  - m2c"
fi
if [[ "$BUILD_TOOLS" -eq "ON" ]]; then
  echo "  - gmsh2aeros"
  echo "  - aeros2dakota"
fi
echo -e "${GREEN}========================================${NO_COLOR}"
echo ""

# -----------------------------------------------------------------------------
# Exit
# -----------------------------------------------------------------------------
exit 0
