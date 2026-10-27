"""Several consensus sequences: the sequences clustered by the MSA columns where they disagree (abPOA's -d).

A column is heterozygous when its second most common letter (a gap counting as one) is carried by at least
`min_freq` of the sequences. Over those columns the sequences are clustered by k-modes (Hamming distance, each
cluster's mode its commonest letter per column), seeded with the sequence nearest all others and then each time the
one furthest from the seeds so far, its distance weighted by how many sequences share its letters there (a lone
read's error is never a seed over an allele many carry); a cluster left smaller than `min_freq` of the sequences
joins its nearest large neighbour.
"""


comptime MAX_ROUNDS = 50
"""k-modes rounds at most (they stop when no row moves)."""


def _hamming(a: List[UInt8], b: List[UInt8]) -> Int:
    var d = 0
    for c in range(len(a)):
        d += Int(a[c] != b[c])
    return d


def _commonest(letters: List[List[UInt8]], label: List[Int], k: Int, mode: List[UInt8]) -> List[UInt8]:
    """Per column the commonest letter of the rows labelled `k` (the mode's own on a tie; `mode` with no rows)."""
    var out = mode.copy()
    var counts = List[Int](length=256, fill=0)
    for c in range(len(mode)):
        for b in range(256):
            counts[b] = 0
        for i in range(len(letters)):
            if label[i] == k:
                counts[Int(letters[i][c])] += 1
        var best = Int(mode[c])
        for b in range(256):
            if counts[b] > counts[best]:
                best = b
        out[c] = UInt8(best)
    return out^


def cluster_rows(rows: List[String], max_clusters: Int, min_freq: Float64) raises -> List[List[Int]]:
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
    var counts = List[Int](length=256, fill=0)
    for c in range(width):
        for b in range(256):
            counts[b] = 0
        for r in rows:
            counts[Int(r.as_bytes()[c])] += 1
        var first = 0
        var second = 0
        for b in range(256):
            if counts[b] > first:
                second = first
                first = counts[b]
            elif counts[b] > second:
                second = counts[b]
        if second >= floor:
            columns.append(c)
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
    for i in range(n):
        var total = 0
        for j in range(n):
            var d = _hamming(letters[i], letters[j])
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
                near = min(near, _hamming(letters[i], m))
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
                var d = _hamming(letters[i], modes[k])
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
                var d = _hamming(letters[i], modes[k])
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
