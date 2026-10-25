"""An MSA anchored on a reference: `on_reference`, each other sequence aligned to the first (CHM13's allele in a
pangenome region) exactly, or pinned as minimap2 pins it (`pins.pinned_pairing`) with the reference's flanks, and its
letters put in the reference letters' columns, so each sequence's alignment with the reference is its own, as a
caller reporting variants against the reference needs."""

from std.atomic import Atomic

from karan import Aligner, Costs, Mode, Ties
from karan.common import spread, thread_count

from .graph import DELETION, INSERTION, cigar_runs
from .pins import pinned_pairing

comptime GAP = UInt8(ord("-"))


comptime REFERENCE_MEMORY = 2 << 30
"""`on_reference`'s memory cap per pair (karan's default, 80 MB, refuses some divergent pairs of a few kb)."""


def on_reference(
    seqs: List[List[UInt8]],
    costs: Costs,
    threads: Int = 1,
    fallback: Optional[Costs] = None,
    flanks: List[List[UInt8]] = List[List[UInt8]](),
) raises -> List[List[UInt8]]:
    """The MSA of `seqs` anchored on the first, the reference: every other distinct sequence aligned to it exactly
    (karan, `costs`, gaps shifted left through repeats as variant callers normalize them) and its letters put in the
    reference letters' columns, the letters it inserts in the columns between them (the runs different sequences
    insert at one place aligned to each other the same way, on the longest). So the alignment of the reference with
    every sequence is that pair's own optimal one: what a caller reports against the reference path keeps each event
    whole, where a placement by many partners' votes can split a copy-number change into copy-sized pieces. Rows in
    `seqs`' order.

    The alignments to the reference run over `threads` threads, an aligner each, the longest sequences first; the
    rows are the same whatever the number.

    With `flanks`, the reference's left and right flanks, each alignment to the reference is pinned as minimap2
    pins it (`pinned_pairing`: its seeds and chain over flank + sequence + flank, exact between the pins), so its
    differences are written as minimap2, and the truth sets dipcall makes with it, write them; no `fallback` (a
    stretch minimap2 deletes and inserts again is its writing too).

    With `fallback`, where at least `REPLACED_SEQUENCES` of the sequences' alignments at `costs` replace a stretch of
    the reference (`_replaces`: a deletion and an insertion of `REPLACEMENT` letters or more side by side) and those
    stretches insert at least `REPLACED_LETTERS` of their letters, those sequences are aligned at `fallback` instead.
    Costs made for close sequences (minimap2's asm5) delete and insert again what lies past their divergence: reads
    of one sequence a few percent apart each are past it nearly everywhere, and costs with dearer long gaps align them
    letter by letter, where the alleles of a population past it at one stretch keep the costs' own records."""
    var n = len(seqs)
    if n == 0:
        return List[List[UInt8]]()
    var distinct = _distinct(seqs)
    var pinned = len(flanks) == 2
    var partner = _pairings(seqs, distinct, costs, threads, flanks)
    # Pinned alignments keep minimap2's writing, a stretch deleted and inserted again included: no fallback.
    if fallback and not pinned:
        _fall_back(seqs, distinct, partner, fallback.value())
    var rows_of = _on_reference_columns(seqs, distinct, partner, costs, fallback)
    var out = List[List[UInt8]]()
    for k in range(n):
        out.append(rows_of[_index_of(seqs, distinct, k)].copy())
    return out^


def _pairings(
    seqs: List[List[UInt8]], distinct: List[Int], costs: Costs, threads: Int, flanks: List[List[UInt8]]
) raises -> List[List[Int]]:
    """For each distinct sequence, the reference letter each of its letters pairs with (-1: inserted): the reference
    itself letter by letter, every other aligned to it (`_onto`, or `pinned_pairing` with two `flanks`). The
    alignments are independent: `threads` workers, an aligner each, take the next, the longest first."""
    var partner = List[List[Int]](length=len(distinct), fill=List[Int]())
    for i in range(len(seqs[distinct[0]])):
        partner[0].append(i)
    var order = List[UInt64]()
    for d in range(1, len(distinct)):
        order.append((UInt64(len(seqs[distinct[d]])) << 32) | UInt64(d))
    sort(order)
    order.reverse()
    var pinned = len(flanks) == 2
    var out = partner.unsafe_ptr()
    var failed = List[Bool](length=len(distinct), fill=False)
    var flags = failed.unsafe_ptr()
    var taken = Atomic[Int64](0)
    var workers = thread_count(threads, len(order))

    def pair(
        slot: Int,
    ) {mut taken, imm order, imm seqs, imm distinct, imm flanks, imm costs, imm pinned, imm out, imm flags}:
        var aligner = Aligner()
        while True:
            var k = Int(taken.fetch_add(1))
            if k >= len(order):
                return
            var d = Int(order[k] & 0xFFFFFFFF)
            try:
                out[unsafe_offset=d] = _pairing_of(aligner, seqs[distinct[0]], seqs[distinct[d]], costs, flanks, pinned)
            except:
                flags[unsafe_offset=d] = True

    spread(pair, workers, workers)
    # One that raised raises here, as one thread would have.
    var aligner = Aligner()
    for d in range(1, len(distinct)):
        if failed[d]:
            partner[d] = _pairing_of(aligner, seqs[distinct[0]], seqs[distinct[d]], costs, flanks, pinned)
    return partner^


def _pairing_of(
    mut aligner: Aligner,
    reference: List[UInt8],
    q: List[UInt8],
    costs: Costs,
    flanks: List[List[UInt8]],
    pinned: Bool,
) raises -> List[Int]:
    """`q`'s pairing with `reference`: pinned between `flanks` when `pinned`, else exact."""
    if pinned:
        return pinned_pairing(aligner, reference, q, flanks[0], flanks[1], costs)
    return _onto(aligner, reference, q, costs)


def _fall_back(seqs: List[List[UInt8]], distinct: List[Int], mut partner: List[List[Int]], fallback: Costs) raises:
    """`on_reference`'s fallback: when most sequences replace stretches of the reference and those stretches are
    much of their letters (reads of one sequence, not alleles; alleles past asm5's divergence replace one stretch of
    a few of them, and their records are written so), those sequences aligned again at `fallback`."""
    var replaced = List[Int]()
    var inserted = 0
    var letters = 0
    for d in range(1, len(distinct)):
        letters += len(seqs[distinct[d]])
        if _replaces(partner[d], len(seqs[distinct[0]])):
            replaced.append(d)
            inserted += _long_inserted(partner[d])
    var most = Float64(len(replaced)) >= REPLACED_SEQUENCES * Float64(len(distinct) - 1)
    if not most or Float64(inserted) < REPLACED_LETTERS * Float64(letters):
        return
    var aligner = Aligner()
    for d in replaced:
        partner[d] = _onto(aligner, seqs[distinct[0]], seqs[distinct[d]], fallback)


comptime REPLACED_SEQUENCES = 0.5
"""The least share of the sequences whose alignment at `on_reference`'s costs replaces a stretch, for its fallback."""
comptime REPLACED_LETTERS = 0.15
"""The least share of the sequences' letters those replacements insert, for `on_reference`'s fallback."""
comptime REPLACEMENT = 50
"""The least length of the deletion and the insertion, within as many letters of each other, by which `_replaces`
finds an alignment that replaced a stretch."""


def _replaces(p: List[Int], reference_length: Int) -> Bool:
    """Whether an alignment (`p` the query letters' reference letters, -1 inserted) deletes `REPLACEMENT` reference
    letters or more and inserts as many query letters within `REPLACEMENT` letters of them: a stretch replaced rather
    than aligned."""
    var inserted_at = List[Int]()
    var deleted_from = List[Int]()
    var deleted_to = List[Int]()
    var next_ref = 0
    var inserting = 0
    for j in range(len(p)):
        if p[j] < 0:
            inserting += 1
            continue
        if inserting >= REPLACEMENT:
            inserted_at.append(next_ref)
        inserting = 0
        if p[j] - next_ref >= REPLACEMENT:
            deleted_from.append(next_ref)
            deleted_to.append(p[j])
        next_ref = p[j] + 1
    if inserting >= REPLACEMENT:
        inserted_at.append(next_ref)
    if reference_length - next_ref >= REPLACEMENT:
        deleted_from.append(next_ref)
        deleted_to.append(reference_length)
    for at in inserted_at:
        for d in range(len(deleted_from)):
            if at >= deleted_from[d] - REPLACEMENT and at <= deleted_to[d] + REPLACEMENT:
                return True
    return False


def _long_inserted(p: List[Int]) -> Int:
    """The query letters in an alignment's insertions of `REPLACEMENT` letters or more (`p` the query letters'
    reference letters, -1 inserted)."""
    var total = 0
    var inserting = 0
    for j in range(len(p) + 1):
        if j < len(p) and p[j] < 0:
            inserting += 1
            continue
        if inserting >= REPLACEMENT:
            total += inserting
        inserting = 0
    return total


def _distinct(seqs: List[List[UInt8]]) -> List[Int]:
    """The first index of each distinct sequence, in order."""
    var distinct = List[Int]()
    for k in range(len(seqs)):
        var seen = False
        for d in range(len(distinct)):
            if seqs[distinct[d]] == seqs[k]:
                seen = True
                break
        if not seen:
            distinct.append(k)
    return distinct^


def _index_of(seqs: List[List[UInt8]], distinct: List[Int], k: Int) -> Int:
    """Which distinct sequence sequence `k` is."""
    for d in range(len(distinct)):
        if seqs[distinct[d]] == seqs[k]:
            return d
    return -1


def _on_reference_columns(
    seqs: List[List[UInt8]],
    distinct: List[Int],
    partner: List[List[Int]],
    costs: Costs,
    fallback: Optional[Costs] = None,
) raises -> List[List[UInt8]]:
    """The rows of the distinct sequences (`distinct[0]` the reference) whose letters `partner` pairs with reference
    letters (-1: inserted): each letter on its reference letter's column, the runs inserted at one place aligned to
    each other on the longest (`on_reference`, `costs`), so their copies share columns."""
    ref c = seqs[distinct[0]]
    var width = len(c)
    # For an inserted letter, its slot (after reference letter i, -1 for before the first) and its place in that
    # slot's run.
    var slot = List[List[Int]]()
    var place_in = List[List[Int]]()
    var slot_width = List[Int](length=width + 1, fill=0)
    for d in range(len(distinct)):
        var slots = List[Int]()
        var places = List[Int]()
        var last = -1
        var run = 0
        for j in range(len(partner[d])):
            if partner[d][j] >= 0:
                slots.append(-2)
                places.append(0)
                last = partner[d][j]
                run = 0
            else:
                slots.append(last)
                places.append(run)
                run += 1
                slot_width[last + 1] = max(slot_width[last + 1], run)
        slot.append(slots^)
        place_in.append(places^)
    # Where several sequences insert different letters at one place, those runs aligned to each other the same way
    # (on the longest of them), so their copies share columns rather than stand left-justified side by side.
    var start_of = List[Dict[Int, Int]]()
    for d in range(len(distinct)):
        var starts = Dict[Int, Int]()
        for j in range(len(partner[d])):
            if partner[d][j] < 0 and (j == 0 or partner[d][j - 1] >= 0 or slot[d][j - 1] != slot[d][j]):
                starts[slot[d][j] + 1] = j
        start_of.append(starts^)
    for place in range(width + 1):
        if slot_width[place] < 2:
            continue
        var members = List[Int]()
        var pieces = List[List[UInt8]]()
        for d in range(len(distinct)):
            var j0 = start_of[d].get(place)
            if not j0:
                continue
            var j = j0.value()
            var piece = List[UInt8]()
            while j < len(partner[d]) and partner[d][j] < 0 and slot[d][j] + 1 == place:
                piece.append(seqs[distinct[d]][j])
                j += 1
            members.append(d)
            pieces.append(piece^)
        var differ = False
        for k in range(1, len(pieces)):
            if pieces[k] != pieces[0]:
                differ = True
                break
        if not differ:
            continue
        # The longest piece first: the anchor of this place's own alignment.
        var longest = 0
        for k in range(1, len(pieces)):
            if len(pieces[k]) > len(pieces[longest]):
                longest = k
        var ordered = List[List[UInt8]]()
        var at = List[Int]()
        ordered.append(pieces[longest].copy())
        at.append(longest)
        for k in range(len(pieces)):
            if k != longest:
                ordered.append(pieces[k].copy())
                at.append(k)
        var sub = on_reference(ordered, costs, fallback=fallback)
        slot_width[place] = len(sub[0])
        for t in range(len(at)):
            var d = members[at[t]]
            var j = start_of[d][place]
            for c in range(len(sub[t])):
                if sub[t][c] != GAP:
                    place_in[d][j] = c
                    j += 1
    # Columns: slot before the first reference letter, then each reference letter followed by its slot.
    var column_of_letter = List[Int](length=width, fill=0)
    var slot_start = List[Int](length=width + 1, fill=0)
    var col = 0
    slot_start[0] = 0
    col += slot_width[0]
    for i in range(width):
        column_of_letter[i] = col
        col += 1
        slot_start[i + 1] = col
        col += slot_width[i + 1]
    var rows_of = List[List[UInt8]]()
    for d in range(len(distinct)):
        ref q = seqs[distinct[d]]
        var row = List[UInt8](length=col, fill=GAP)
        for j in range(len(q)):
            if partner[d][j] >= 0:
                row[column_of_letter[partner[d][j]]] = q[j]
            else:
                row[slot_start[slot[d][j] + 1] + place_in[d][j]] = q[j]
        rows_of.append(row^)
    return rows_of^


def _onto(mut aligner: Aligner, c: List[UInt8], q: List[UInt8], costs: Costs) raises -> List[Int]:
    """For each letter of `q`, the letter of `c` its exact alignment to `c` pairs it with (-1: inserted)."""
    var found = aligner.align(Span(c), Span(q), costs, Mode.GLOBAL, ties=Ties.LEFT, max_memory=REFERENCE_MEMORY)
    return _pairing(found.cigar, len(q))


def _pairing(cigar: String, q_length: Int) raises -> List[Int]:
    """For each query letter of a global alignment's CIGAR, the reference letter it pairs with (-1: inserted)."""
    var runs = cigar_runs(cigar)
    var out = List[Int](capacity=q_length)
    var i = 0
    for x in range(len(runs.ops)):
        var op = runs.ops[x]
        for _ in range(runs.counts[x]):
            if op == DELETION:
                i += 1
            elif op == INSERTION:
                out.append(-1)
            else:
                out.append(i)
                i += 1
    return out^
