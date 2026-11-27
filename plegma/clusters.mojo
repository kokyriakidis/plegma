"""Several consensus sequences: the sequences clustered by the MSA columns where they disagree (abPOA's -d).

A column is heterozygous when its second most common letter (a gap counting as one) is carried by at least
`min_freq` of the sequences. Over those columns the sequences are clustered by k-modes (Hamming distance, a column of
two letters counting `SUBSTITUTION_WEIGHT` columns where a gap is one of the two; each cluster's mode its commonest
letter per column), seeded with the sequence nearest all others and then each time the
one furthest from the seeds so far, its distance weighted by how many sequences share its letters there (a lone
read's error is never a seed over an allele many carry); a cluster left smaller than `min_freq` of the sequences
joins its nearest large neighbour. A row's `MISSING` letters, where it has no letter to give (a read anchored at one
end, past its far end), count for no column and set it apart from none.
"""


from std.atomic import Atomic

from karan.common import spread, thread_count


comptime MAX_ROUNDS = 50
"""k-modes rounds at most (they stop when no row moves)."""

comptime GAP = UInt8(ord("-"))

comptime MISSING = UInt8(ord("."))
"""A row's letter where it covers no column (`cluster_rows`)."""

comptime SUBSTITUTION_WEIGHT = 3
"""A column where the rows disagree by two letters counts this many times one where some rows have a gap: HiFi reads
slip in runs (a gap) far more often than they misread a letter, so the letters apart are the likelier alleles."""


def _hamming(a: List[UInt8], b: List[UInt8], weights: List[Int]) -> Int:
    var d = 0
    for c in range(len(a)):
        if a[c] != b[c] and a[c] != MISSING and b[c] != MISSING:
            d += weights[c]
    return d


def _commonest(letters: List[List[UInt8]], label: List[Int], k: Int, mode: List[UInt8]) -> List[UInt8]:
    """Per column the commonest letter of the rows labelled `k` (the mode's own on a tie; `mode` with no rows)."""
    var out = mode.copy()
    var counts = List[Int](length=256, fill=0)
    var touched = List[Int](capacity=8)
    for c in range(len(mode)):
        touched.clear()
        for i in range(len(letters)):
            if label[i] == k:
                var b = Int(letters[i][c])
                if b == Int(MISSING):
                    continue
                if counts[b] == 0:
                    touched.append(b)
                counts[b] += 1
        # Keep ascending byte order on ties (identical to scanning `b in range(256)`).
        for a in range(1, len(touched)):
            var key = touched[a]
            var j = a
            while j > 0 and touched[j - 1] > key:
                touched[j] = touched[j - 1]
                j -= 1
            touched[j] = key
        var best = Int(mode[c])
        for b in touched:
            if counts[b] > counts[best]:
                best = b
        for b in touched:
            counts[b] = 0
        out[c] = UInt8(best)
    return out^


def cluster_rows(rows: List[String], max_clusters: Int, min_freq: Float64, threads: Int = 1) raises -> List[List[Int]]:
    """The rows' clusters, at most `max_clusters`, each as row indices, in order of their first rows."""
    var n = len(rows)
    var everyone = List[Int]()
    for i in range(n):
        everyone.append(i)
    var out = List[List[Int]]()
    if n == 0 or max_clusters < 2:
        out.append(everyone^)
        return out^
    var width = rows[0].byte_length()
    var floor = max(1, Int(min_freq * Float64(n) + 0.999999))
    var columns = List[Int]()
    var weights = List[Int]()
    var counts = List[Int](length=256, fill=0)
    var touched = List[Int](capacity=8)
    for c in range(width):
        touched.clear()
        for r in rows:
            var b = Int(r.as_bytes()[c])
            if b == Int(MISSING):
                continue
            if counts[b] == 0:
                touched.append(b)
            counts[b] += 1
        if len(touched) >= 2:
            for a in range(1, len(touched)):
                var key = touched[a]
                var j = a
                while j > 0 and touched[j - 1] > key:
                    touched[j] = touched[j - 1]
                    j -= 1
                touched[j] = key
            var first = 0
            var second = 0
            var first_letter = 0
            var second_letter = 0
            for b in touched:
                if counts[b] > first:
                    second = first
                    second_letter = first_letter
                    first = counts[b]
                    first_letter = b
                elif counts[b] > second:
                    second = counts[b]
                    second_letter = b
            if second >= floor:
                columns.append(c)
                var substitution = UInt8(first_letter) != GAP and UInt8(second_letter) != GAP
                weights.append(SUBSTITUTION_WEIGHT if substitution else 1)
        for b in touched:
            counts[b] = 0
    if len(columns) == 0:
        out.append(everyone^)
        return out^
    # k-modes over those columns, which a few reads' errors do not move: the first mode the row nearest every other
    # (least summed Hamming distance), each next the row furthest from the modes so far times the rows sharing its
    # letters; then every row joins its nearest mode and each mode becomes its rows' commonest letter per column,
    # until no row moves. Unweighted, one read's third letter in the only differing column (chr18:25,307,645, 38
    # reads A, 38 a gap, one G) seeded the second mode, every gap row tied between the two and joined the first, and
    # the het was lost.
    var letters = List[List[UInt8]](capacity=n)
    for i in range(n):
        var key = List[UInt8](capacity=len(columns))
        for c in columns:
            key.append(rows[i].as_bytes()[c])
        letters.append(key^)
    var first_mode = 0
    var least = Int.MAX
    var sharing = List[Int](length=n, fill=0)
    var workers = thread_count(threads, n)
    if workers > 1 and n >= 32 and len(columns) >= 4:
        var totals = List[Int](length=n, fill=0)
        var totals_p = totals.unsafe_ptr()
        var sharing_p = sharing.unsafe_ptr()
        var taken = Atomic[Int64](0)

        def row_distances(slot: Int) {mut taken, imm n, imm letters, imm weights, imm totals_p, imm sharing_p}:
            while True:
                var i = Int(taken.fetch_add(1))
                if i >= n:
                    return
                var total = 0
                var shared = 0
                for j in range(n):
                    var d = _hamming(letters[i], letters[j], weights)
                    total += d
                    shared += Int(d == 0)
                totals_p[unsafe_offset=i] = total
                sharing_p[unsafe_offset=i] = shared

        spread(row_distances, workers, workers)
        for i in range(n):
            if totals[i] < least:
                least = totals[i]
                first_mode = i
    else:
        for i in range(n):
            var total = 0
            for j in range(n):
                var d = _hamming(letters[i], letters[j], weights)
                total += d
                sharing[i] += Int(d == 0)
            if total < least:
                least = total
                first_mode = i
    var modes = List[List[UInt8]]()
    modes.append(letters[first_mode].copy())
    while len(modes) < max_clusters:
        var furthest = -1
        var far = 0
        for i in range(n):
            var near = Int.MAX
            for m in modes:
                near = min(near, _hamming(letters[i], m, weights))
            if near * sharing[i] > far:
                far = near * sharing[i]
                furthest = i
        if furthest < 0:
            break
        modes.append(letters[furthest].copy())
    if len(modes) < 2:
        out.append(everyone^)
        return out^
    var label = List[Int](length=n, fill=-1)
    for _ in range(MAX_ROUNDS):
        var moved = False
        for i in range(n):
            var best = 0
            var best_d = Int.MAX
            for k in range(len(modes)):
                var d = _hamming(letters[i], modes[k], weights)
                if d < best_d:
                    best_d = d
                    best = k
            if best != label[i]:
                moved = True
                label[i] = best
        if not moved:
            break
        for k in range(len(modes)):
            modes[k] = _commonest(letters, label, k, modes[k])
    # Clusters under the floor join the nearest large cluster's mode.
    var size = List[Int](length=len(modes), fill=0)
    for i in range(n):
        size[label[i]] += 1
    var large = 0
    for k in range(len(modes)):
        if size[k] >= floor:
            large += 1
    if large < 2:
        out.append(everyone^)
        return out^
    for i in range(n):
        if size[label[i]] >= floor:
            continue
        var best = label[i]
        var best_d = Int.MAX
        for k in range(len(modes)):
            if size[k] >= floor:
                var d = _hamming(letters[i], modes[k], weights)
                if d < best_d:
                    best_d = d
                    best = k
        label[i] = best
    for k in range(len(modes)):
        var members = List[Int]()
        for i in range(n):
            if label[i] == k:
                members.append(i)
        if len(members) > 0:
            out.append(members^)
    # In order of their first rows, as abPOA lists them.
    for a in range(1, len(out)):
        var at = a
        while at > 0 and out[at - 1][0] > out[at][0]:
            out.swap_elements(at - 1, at)
            at -= 1
    return out^
