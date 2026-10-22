"""Several consensus sequences: the sequences clustered by the MSA columns where they disagree (abPOA's -d).

A column is heterozygous when its second most common letter (a gap counting as one) is carried by at least
`min_freq` of the sequences. The medoids are the most common patterns of letters over those columns that at
least `min_freq` of the sequences carry, and every sequence joins the nearest at Hamming distance; a cluster
left smaller than `min_freq` of the sequences joins its nearest large neighbour.
"""

from std.collections import Dict


def _distance(rows: List[String], a: Int, b: Int, columns: List[Int]) -> Int:
    var d = 0
    var x = rows[a].as_bytes()
    var y = rows[b].as_bytes()
    for c in columns:
        d += Int(x[c] != y[c])
    return d


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
    # Medoids: the most common patterns of letters over those columns, each carried by at least `floor` rows
    # (noise in a few rows cannot make a medoid), most common first; every row goes to its nearest one.
    var patterns = Dict[String, Int]()
    var first_row = Dict[String, Int]()
    for i in range(n):
        var key = List[UInt8](capacity=len(columns))
        for c in columns:
            key.append(rows[i].as_bytes()[c])
        var text = String(unsafe_from_utf8=key)
        if text in patterns:
            patterns[text] += 1
        else:
            patterns[text] = 1
            first_row[text] = i
    var candidates = List[Tuple[Int, Int]]()  # (count, row)
    for entry in patterns.items():
        if entry.value >= floor:
            candidates.append((entry.value, first_row[entry.key]))
    var medoids = List[Int]()
    while len(medoids) < max_clusters and len(candidates) > 0:
        var best = 0
        for k in range(1, len(candidates)):
            if candidates[k][0] > candidates[best][0] or (
                candidates[k][0] == candidates[best][0] and candidates[k][1] < candidates[best][1]
            ):
                best = k
        medoids.append(candidates[best][1])
        _ = candidates.pop(best)
    if len(medoids) < 2:
        out.append(everyone^)
        return out^
    var label = List[Int](length=n, fill=0)
    for i in range(n):
        var best = 0
        var best_d = Int.MAX
        for k in range(len(medoids)):
            var d = _distance(rows, i, medoids[k], columns)
            if d < best_d:
                best_d = d
                best = k
        label[i] = best
    # Clusters under the floor join the nearest large cluster's medoid.
    var size = List[Int](length=len(medoids), fill=0)
    for i in range(n):
        size[label[i]] += 1
    for i in range(n):
        if size[label[i]] >= floor:
            continue
        var best = label[i]
        var best_d = Int.MAX
        for k in range(len(medoids)):
            if size[k] >= floor:
                var d = _distance(rows, i, medoids[k], columns)
                if d < best_d:
                    best_d = d
                    best = k
        label[i] = best
    for k in range(len(medoids)):
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
