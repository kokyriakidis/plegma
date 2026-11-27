"""An MSA anchored on a reference: `on_reference`, each other sequence aligned to the first (CHM13's allele in a
pangenome region) exactly, or pinned as minimap2 pins it (`pins.pinned_pairing_with_index`) with the reference's
flanks, and its letters put in the reference letters' columns, so each sequence's alignment with the reference is
its own, as a caller reporting variants against the reference needs."""

from std.atomic import Atomic

from karan import Aligner, Costs, Mode, Ties
from karan.common import spread, thread_count

from .cigar import pairing_from_cigar
from .pins import ReferenceIndex, pinned_pairing_with_index
from .seeds import POSITION_BITS
from .trace import REFERENCE, Trace

comptime GAP = UInt8(ord("-"))
comptime REFERENCE_MEMORY = 2 << 30
"""`on_reference`'s memory cap per pair (karan's default, 80 MB, refuses some divergent pairs of a few kb)."""
comptime REPLACED_SEQUENCES = 0.5
"""The least share of the sequences whose alignment at `on_reference`'s costs replaces a stretch, for its fallback."""
comptime REPLACED_LETTERS = 0.15
"""The least share of the sequences' letters those replacements insert, for `on_reference`'s fallback."""
comptime REPLACEMENT = 50
"""The least length of the deletion and the insertion, within as many letters of each other, by which `_replaces`
finds an alignment that replaced a stretch."""


def _concat3(a: List[UInt8], b: List[UInt8], c: List[UInt8]) -> List[UInt8]:
    var out = List[UInt8](capacity=len(a) + len(b) + len(c))
    for x in a:
        out.append(x)
    for x in b:
        out.append(x)
    for x in c:
        out.append(x)
    return out^


def _dedup(seqs: List[List[UInt8]]) -> Tuple[List[Int], List[Int]]:
    """Deduplicates `seqs` in O(N) time: returns `(distinct, index_of)` where `distinct` holds the first index of
    each unique sequence in order and `index_of[k]` maps sequence `k` to its slot in `distinct`."""
    var distinct = List[Int]()
    var index_of = List[Int](capacity=len(seqs))
    var seen = Dict[String, Int]()
    for k in range(len(seqs)):
        var key = String(unsafe_from_utf8=seqs[k])
        var d = seen.get(key, -1)
        if d < 0:
            d = len(distinct)
            seen[key] = d
            distinct.append(k)
        index_of.append(d)
    return (distinct^, index_of^)


def on_reference(
    seqs: List[List[UInt8]],
    costs: Costs,
    threads: Int = 1,
    fallback: Optional[Costs] = None,
    flanks: List[List[UInt8]] = List[List[UInt8]](),
) raises -> List[List[UInt8]]:
    """The MSA of `seqs` anchored on the first, the reference: every other distinct sequence aligned to it exactly
    (or pinned with `flanks`) in parallel over `threads`, fallback re-aligned in parallel when triggered, and
    multi-sequence insertion slots aligned in parallel."""
    var n = len(seqs)
    if n == 0:
        return List[List[UInt8]]()
    var deduped = _dedup(seqs)
    ref distinct = deduped[0]
    ref index_of = deduped[1]
    var pinned = len(flanks) == 2
    var trace = Trace()
    if trace.on(REFERENCE):
        trace.line(
            REFERENCE,
            String(
                n,
                " sequences (",
                len(distinct),
                " distinct), reference ",
                len(seqs[0]),
                " letters; ",
                String(
                    "pinned between flanks of ", len(flanks[0]), " and ", len(flanks[1]), " letters"
                ) if pinned else String("exact"),
                ", ",
                threads,
                " threads",
            ),
        )
    var partner = _pairings(seqs, distinct, costs, threads, flanks)
    if fallback and not pinned:
        _fall_back(seqs, distinct, partner, fallback.value(), threads)
    var rows_of = _on_reference_columns(seqs, distinct, partner, costs, None if pinned else fallback, threads)
    var out = List[List[UInt8]](capacity=n)
    for k in range(n):
        out.append(rows_of[index_of[k]].copy())
    return out^


def _pairings(
    seqs: List[List[UInt8]], distinct: List[Int], costs: Costs, threads: Int, flanks: List[List[UInt8]]
) raises -> List[List[Int]]:
    """For each distinct sequence, the reference letter each of its letters pairs with (-1: inserted). When pinned,
    builds the reference minimizer index once and shares it read-only across all worker threads."""
    var partner = List[List[Int]](length=len(distinct), fill=List[Int]())
    ref reference = seqs[distinct[0]]
    for i in range(len(reference)):
        partner[0].append(i)
    if len(distinct) <= 1:
        return partner^
    var order = List[UInt64](capacity=len(distinct) - 1)
    for d in range(1, len(distinct)):
        order.append((UInt64(len(seqs[distinct[d]])) << 32) | UInt64(d))
    sort(order)
    order.reverse()
    var pinned = len(flanks) == 2
    if pinned:
        # Check reach up front before indexing so oversized sequences fail immediately with the exact error.
        for d in range(1, len(distinct)):
            var longest = len(flanks[0]) + max(len(reference), len(seqs[distinct[d]])) + len(flanks[1])
            if longest >= 1 << POSITION_BITS:
                raise Error(
                    t"pinned alignment: {longest} letters with the flanks, past the 2^{POSITION_BITS} pins reach"
                )
        var whole_ref = _concat3(flanks[0], reference, flanks[1])
        var ref_index = ReferenceIndex(whole_ref)
        var workers = thread_count(threads, len(order))
        if workers <= 1:
            # Single sequence or single thread: let `pinned_pairing_with_index` parallelize its windows directly.
            var aligner = Aligner(warm=False)
            for k in range(len(order)):
                var d = Int(order[k] & 0xFFFFFFFF)
                partner[d] = pinned_pairing_with_index(
                    aligner, ref_index, reference, seqs[distinct[d]], flanks[0], flanks[1], costs, threads
                )
            return partner^
        # Multiple distinct sequences: keep inner_threads = 1 to avoid nested `spread` oversubscription/deadlock.
        var out = partner.unsafe_ptr()
        var failed = List[Bool](length=len(distinct), fill=False)
        var flags = failed.unsafe_ptr()
        var taken = Atomic[Int64](0)

        def pair_pinned(
            slot: Int,
        ) {mut taken, imm order, imm seqs, imm distinct, imm ref_index, imm flanks, imm costs, imm out, imm flags,}:
            var aligner = Aligner(warm=False)
            while True:
                var k = Int(taken.fetch_add(1))
                if k >= len(order):
                    return
                var d = Int(order[k] & 0xFFFFFFFF)
                try:
                    out[unsafe_offset=d] = pinned_pairing_with_index(
                        aligner,
                        ref_index,
                        seqs[distinct[0]],
                        seqs[distinct[d]],
                        flanks[0],
                        flanks[1],
                        costs,
                        1,
                    )
                except:
                    flags[unsafe_offset=d] = True

        spread(pair_pinned, workers, workers)
        var any_failed = False
        for d in range(1, len(distinct)):
            if failed[d]:
                any_failed = True
                break
        if any_failed:
            var aligner = Aligner(warm=False)
            for d in range(1, len(distinct)):
                if failed[d]:
                    partner[d] = pinned_pairing_with_index(
                        aligner, ref_index, reference, seqs[distinct[d]], flanks[0], flanks[1], costs
                    )
        return partner^

    var workers = thread_count(threads, len(order))
    var out = partner.unsafe_ptr()
    var failed = List[Bool](length=len(distinct), fill=False)
    var flags = failed.unsafe_ptr()
    var taken = Atomic[Int64](0)

    def pair_exact(
        slot: Int,
    ) {mut taken, imm order, imm seqs, imm distinct, imm costs, imm out, imm flags}:
        var aligner = Aligner(warm=False)
        while True:
            var k = Int(taken.fetch_add(1))
            if k >= len(order):
                return
            var d = Int(order[k] & 0xFFFFFFFF)
            try:
                out[unsafe_offset=d] = _onto(aligner, seqs[distinct[0]], seqs[distinct[d]], costs)
            except:
                flags[unsafe_offset=d] = True

    spread(pair_exact, workers, workers)
    var any_failed = False
    for d in range(1, len(distinct)):
        if failed[d]:
            any_failed = True
            break
    if any_failed:
        var aligner = Aligner(warm=False)
        for d in range(1, len(distinct)):
            if failed[d]:
                partner[d] = _onto(aligner, reference, seqs[distinct[d]], costs)
    return partner^


def _fall_back(
    seqs: List[List[UInt8]], distinct: List[Int], mut partner: List[List[Int]], fallback: Costs, threads: Int = 1
) raises:
    """`on_reference`'s fallback, re-aligning replaced sequences in parallel over `threads`."""
    var replaced = List[Int]()
    var inserted = 0
    var letters = 0
    for d in range(1, len(distinct)):
        letters += len(seqs[distinct[d]])
        if _replaces(partner[d], len(seqs[distinct[0]])):
            replaced.append(d)
            inserted += _long_inserted(partner[d])
    if len(distinct) <= 1:
        return
    var most = Float64(len(replaced)) >= REPLACED_SEQUENCES * Float64(len(distinct) - 1)
    var trace = Trace()
    if not most or Float64(inserted) < REPLACED_LETTERS * Float64(letters):
        if trace.on(REFERENCE) and len(replaced) > 0:
            trace.line(
                REFERENCE,
                String(
                    "no fallback: ",
                    len(replaced),
                    " of ",
                    len(distinct) - 1,
                    " replace a stretch, inserting ",
                    inserted,
                    " of ",
                    letters,
                    " letters",
                ),
            )
        return
    if trace.on(REFERENCE):
        trace.line(
            REFERENCE,
            String(
                "fallback costs for ",
                len(replaced),
                " of ",
                len(distinct) - 1,
                " (stretches replaced, ",
                inserted,
                " of ",
                letters,
                " letters inserted)",
            ),
        )
    var workers = thread_count(threads, len(replaced))
    if workers > 1 and len(replaced) >= 2:
        var order = List[UInt64](capacity=len(replaced))
        for d in replaced:
            order.append((UInt64(len(seqs[distinct[d]])) << 32) | UInt64(d))
        sort(order)
        order.reverse()
        var out = partner.unsafe_ptr()
        var failed = List[Bool](length=len(distinct), fill=False)
        var flags = failed.unsafe_ptr()
        var taken = Atomic[Int64](0)

        def fallback_worker(
            slot: Int,
        ) {mut taken, imm order, imm seqs, imm distinct, imm fallback, imm out, imm flags}:
            var aligner = Aligner(warm=False)
            while True:
                var k = Int(taken.fetch_add(1))
                if k >= len(order):
                    return
                var d = Int(order[k] & 0xFFFFFFFF)
                try:
                    out[unsafe_offset=d] = _onto(aligner, seqs[distinct[0]], seqs[distinct[d]], fallback)
                except:
                    flags[unsafe_offset=d] = True

        spread(fallback_worker, workers, workers)
        var any_failed = False
        for d in replaced:
            if failed[d]:
                any_failed = True
                break
        if any_failed:
            var aligner = Aligner(warm=False)
            for d in replaced:
                if failed[d]:
                    partner[d] = _onto(aligner, seqs[distinct[0]], seqs[distinct[d]], fallback)
        return

    var aligner = Aligner(warm=False)
    for d in replaced:
        partner[d] = _onto(aligner, seqs[distinct[0]], seqs[distinct[d]], fallback)


def _replaces(p: List[Int], reference_length: Int) -> Bool:
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


struct _SlotJob(Copyable, Movable):
    """One insertion slot where multiple sequences insert differing runs that must be sub-aligned."""

    var place: Int
    var members: List[Int]
    var ordered: List[List[UInt8]]
    var at: List[Int]

    def __init__(out self, place: Int, var members: List[Int], var ordered: List[List[UInt8]], var at: List[Int]):
        self.place = place
        self.members = members^
        self.ordered = ordered^
        self.at = at^


def _on_reference_columns(
    seqs: List[List[UInt8]],
    distinct: List[Int],
    partner: List[List[Int]],
    costs: Costs,
    fallback: Optional[Costs] = None,
    threads: Int = 1,
) raises -> List[List[UInt8]]:
    """The rows of the distinct sequences whose letters `partner` pairs with reference letters (-1: inserted):
    multi-sequence differing insertion slots are sub-aligned in parallel over `threads`."""
    ref c = seqs[distinct[0]]
    var width = len(c)
    var slot = List[List[Int]](capacity=len(distinct))
    var place_in = List[List[Int]](capacity=len(distinct))
    var slot_width = List[Int](length=width + 1, fill=0)
    for d in range(len(distinct)):
        var slots = List[Int](capacity=len(partner[d]))
        var places = List[Int](capacity=len(partner[d]))
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
    var start_of = List[Dict[Int, Int]](capacity=len(distinct))
    for d in range(len(distinct)):
        var starts = Dict[Int, Int]()
        for j in range(len(partner[d])):
            if partner[d][j] < 0 and (j == 0 or partner[d][j - 1] >= 0 or slot[d][j - 1] != slot[d][j]):
                starts[slot[d][j] + 1] = j
        start_of.append(starts^)

    # Collect all independent insertion-slot sub-alignment jobs.
    var jobs = List[_SlotJob]()
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
        var longest = 0
        for k in range(1, len(pieces)):
            if len(pieces[k]) > len(pieces[longest]):
                longest = k
        var ordered = List[List[UInt8]](capacity=len(pieces))
        var at = List[Int](capacity=len(pieces))
        ordered.append(pieces[longest].copy())
        at.append(longest)
        for k in range(len(pieces)):
            if k != longest:
                ordered.append(pieces[k].copy())
                at.append(k)
        jobs.append(_SlotJob(place, members^, ordered^, at^))

    var n_jobs = len(jobs)
    var slot_workers = thread_count(threads, n_jobs)
    if slot_workers > 1 and n_jobs >= 2:
        var sub_threads = 1
        var subs = List[List[List[UInt8]]](length=n_jobs, fill=List[List[UInt8]]())
        var failed = List[Bool](length=n_jobs, fill=False)
        var subs_p = subs.unsafe_ptr()
        var failed_p = failed.unsafe_ptr()
        var taken = Atomic[Int64](0)

        def run_slot(
            slot_i: Int,
        ) {mut taken, imm n_jobs, imm jobs, imm costs, imm fallback, imm sub_threads, imm subs_p, imm failed_p}:
            while True:
                var idx = Int(taken.fetch_add(1))
                if idx >= n_jobs:
                    return
                try:
                    subs_p[unsafe_offset=idx] = on_reference(jobs[idx].ordered, costs, sub_threads, fallback=fallback)
                except:
                    failed_p[unsafe_offset=idx] = True

        spread(run_slot, slot_workers, slot_workers)
        for idx in range(n_jobs):
            var sub = on_reference(jobs[idx].ordered, costs, threads, fallback=fallback) if failed[idx] else subs[
                idx
            ].copy()
            ref job = jobs[idx]
            slot_width[job.place] = len(sub[0])
            for t in range(len(job.at)):
                var d = job.members[job.at[t]]
                var j = start_of[d][job.place]
                for col_i in range(len(sub[t])):
                    if sub[t][col_i] != GAP:
                        place_in[d][j] = col_i
                        j += 1
    else:
        for idx in range(n_jobs):
            ref job = jobs[idx]
            var sub = on_reference(job.ordered, costs, threads, fallback=fallback)
            slot_width[job.place] = len(sub[0])
            for t in range(len(job.at)):
                var d = job.members[job.at[t]]
                var j = start_of[d][job.place]
                for col_i in range(len(sub[t])):
                    if sub[t][col_i] != GAP:
                        place_in[d][j] = col_i
                        j += 1

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
    var rows_of = List[List[UInt8]](capacity=len(distinct))
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
    return pairing_from_cigar(found.cigar, len(q))
