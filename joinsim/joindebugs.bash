#!/bin/bash
OUTFILE=$1
shift # Shifts command line arguments down, removing $1
cat "$@" | sort -n | sed -e "s/^\<[0-9][0-9]* //" > $OUTFILE
#cat "$@" | sort -n > $OUTFILE
