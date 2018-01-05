#!/bin/sh

## Verbose mode
test "x$VERBOSE" = xx && set -x

SRC_PREFIX=$srcdir/SFS_small
OUT_PREFIX=$srcdir/out
TEST_ID=out1
OPTIONS=-n 10
INTERPOP=

../src/selink -h 2>/dev/null || exit 1
