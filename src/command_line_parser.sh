#!/bin/bash

# parse command-line arguments using "getopt"
SHORT_OPTIONS=""
LONG_OPTIONS="help,config:,pre:,post:"

PARSER="$(getopt -o "$SHORT_OPTIONS" -l "$LONG_OPTIONS" -- "$@")" || {
  # getopt failed
  printf "*** Error: Could not understand the "
  printf "command line argument \"%s\".\n" "$?"
  exit 1
}

eval set -- "$PARSER"

# Get user specified options.
while true; do

  case "$1" in

    --help)
      echo "Usage: driver.sh"
      echo "  Options:"
      echo "    --config: "
      echo "    --param: "
      echo "    --result: "
      echo "    --pre: "
      echo "    --post: "
      exit 0
      ;;
    --config)
      USER_CONFIG="$2"
      shift 2
      ;;
    --pre)
      PREPROCESS_FILE="$2"
      shift 2
      ;;
    --post)
      POSTPROCESS_FILE="$2"
      shitf 2
      ;;
    --)
      shift
      break
      ;;
    *)
      shift
      break
      ;;

  esac

done
