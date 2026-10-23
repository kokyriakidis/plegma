#!/bin/sh
# abPOA as vntr-eval's poa_abpoa runs it, for candidates.sh: the sequences longest first (ties in input order),
# then abpoa -m 0 with the remaining flags; the last argument is the FASTA.
#   sh experiments/vntreval/abpoa_longest.sh <abpoa> [flags...] <in.fa>
abpoa=$1
shift
eval "in=\${$#}"
flags=""
while [ $# -gt 1 ]; do flags="$flags $1"; shift; done
sorted=$(mktemp)
python3 -I -c '
import sys
names, seqs = [], []
for line in open(sys.argv[1]):
    line = line.rstrip("\n")
    if line.startswith(">"):
        names.append(line); seqs.append([])
    elif line:
        seqs[-1].append(line)
order = sorted(range(len(names)), key=lambda k: (-sum(len(x) for x in seqs[k]), k))
with open(sys.argv[2], "w") as f:
    for k in order:
        f.write(names[k] + "\n" + "".join(seqs[k]) + "\n")
' "$in" "$sorted"
"$abpoa" -m 0 $flags "$sorted"
status=$?
rm -f "$sorted"
exit $status
