#!/usr/bin/env python3
"""Which features of an anchor tell an off one (pairing letters no optimal alignment pairs) from a sound one.

usage: python3 -I experiments/anchors/analyze.py <anchors.tsv> [plan] [--min-depth d]

Over the anchors of one plan (the first in the file by default) at depth `d` or deeper (1: inside the windows the
top level leaves): the share off (an anchor with any letter off, and its letters) by length, by its string's copies,
by its window's size and length difference, by how many anchors its window chained, by repeat-mask cover, and by
where its diagonal lies between the window's two ends: `diag_start` + `diag_end` is the window's length difference
D, the indel an alignment of the window must make somewhere; an anchor at diag_start 0 puts it all after itself, at
D all before, between splits it.
"""
import collections
import csv
import sys


def bucket(x, edges):
    for e in edges:
        if x <= e:
            return "<=%s" % e
    return ">%s" % edges[-1]


def main():
    args = sys.argv[1:]
    min_depth = 1
    if "--min-depth" in args:
        k = args.index("--min-depth")
        min_depth = int(args[k + 1])
        del args[k:k + 2]
    rows = list(csv.DictReader(open(args[0]), delimiter="\t"))
    plan = args[1] if len(args) > 1 else rows[0]["plan"]
    rows = [r for r in rows if r["plan"] == plan and int(r["depth"]) >= min_depth]
    print("plan %s, depth >= %d: %d anchors" % (plan, min_depth, len(rows)))

    def table(name, key):
        groups = collections.defaultdict(lambda: [0, 0, 0, 0])
        for r in rows:
            g = groups[key(r)]
            g[0] += 1
            g[1] += int(r["off"]) > 0
            g[2] += int(r["letters"])
            g[3] += int(r["off"])
        print("\n%s" % name)
        print("  %-16s %8s %8s %9s %9s %7s" % ("value", "anchors", "off", "letters", "off", "off %"))
        for k in sorted(groups, key=lambda v: (len(str(v)), str(v))):
            g = groups[k]
            print("  %-16s %8d %8d %9d %9d %6.1f%%" % (k, g[0], g[1], g[2], g[3], 100.0 * g[3] / max(1, g[2])))

    table("length", lambda r: bucket(int(r["length"]), [25, 35, 50, 75, 100, 200, 400]))
    table("copies (partner x query)", lambda r: "%sx%s" % (r["ref_copies"], r["query_copies"]))
    table("window (larger side)", lambda r: bucket(max(int(r["win_ref"]), int(r["win_query"])),
                                                    [100, 300, 1000, 3000, 10000]))
    table("window length difference |D|", lambda r: bucket(abs(int(r["win_query"]) - int(r["win_ref"])),
                                                            [0, 5, 30, 100, 300, 1000]))
    table("anchors chained in its window", lambda r: bucket(int(r["siblings"]), [1, 2, 5, 10, 30]))
    table("in a masked repeat", lambda r: r["in_repeat"])

    def split(r):
        d = int(r["diag_start"]) + int(r["diag_end"])
        s = int(r["diag_start"])
        if d == 0:
            return "D=0, at 0" if s == 0 else "D=0, off 0"
        f = s / d
        if abs(s) <= 2:
            return "at start (0)"
        if abs(s - d) <= 2:
            return "at end (D)"
        if 0 < f < 1:
            return "between"
        return "outside"

    table("diagonal between the window's ends", split)
    table("depth", lambda r: r["depth"])
    table("level k", lambda r: r["k"])


if __name__ == "__main__":
    main()
