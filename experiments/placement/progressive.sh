#!/bin/sh
# The placement harness's progressive mode called as iterate.py calls plegma, its MSA on stdout:
#   progressive.sh <strategy> <refine passes> [group band[,max band]] -r 1 <in.fa>
# PLACEMENT names the harness binary (build/placement by default); the group band (0, off) turns on the refinement's
# group moves, the max band (512) caps how far its bands widen.
strategy=$1
passes=$2
shift 2
group=0
if [ "$1" != "-r" ]; then
    group=$1
    shift
fi
in=$3
out=$(mktemp)
${PLACEMENT:-$(dirname "$0")/../../build/placement} --progressive "$in" "$out" "$strategy" "$passes" "$group" || { rm -f "$out"; exit 1; }
cat "$out"
rm -f "$out"
