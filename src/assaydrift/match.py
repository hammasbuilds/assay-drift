"""Does this primer still bind this genome?

The whole study reduces to that question, asked once per oligo per sequence.
Three details decide whether the answer means anything.

**An `N` in the target is missing data, not a mismatch.** Sequencing gaps are
written as `N`, and a genome with an `N` under a primer binding site simply does
not say whether the primer matches. Counting those as mismatches would be the
single worst mistake available here: sequencing quality changed enormously over
the pandemic, so it would manufacture a *time trend* out of laboratory practice
and present it as viral drift. Oligos whose binding site contains any `N` are
reported separately and excluded from the rates, because "unknown" is not
"broken" and it is not "fine" either.

**Ambiguity codes in the primer are deliberate.** `R` means the oligo was
synthesised as a mixture of A and G, so it genuinely matches both. The Charite
RdRp assay uses them on purpose, to catch SARS-related viruses beyond SARS-CoV-2.
Treating `R` as a literal letter would report that assay as broken everywhere.

**A mismatch at the 3' end is not like a mismatch in the middle.** Polymerase
extends from the 3' terminus, so a mismatch in the last few bases can stop
amplification outright while the same mismatch in the middle is often tolerated.
They are counted separately rather than summed into one number that hides the
difference.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace

# What each IUPAC code is allowed to be.
IUPAC: dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}

COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

# How many bases at the 3' end count as "the business end" of a primer.
# Five is the usual rule of thumb in primer design: a mismatch inside this
# window is the kind that stops extension rather than merely slowing it.
THREE_PRIME_WINDOW = 5

# Above this many definite mismatches, an alignment is not the primer's binding
# site any more - it is noise found somewhere else in a 30,000 base genome.
# The absolute ceiling, kept for the longest oligos and for the NOT_FOUND sentinel.
MAX_MISMATCHES = 8

# An absolute cap cannot do this job, because the catalogue's oligos are 18 to 26 bases
# and 8 mismatches is 31% of a 26-mer but 44% of an 18-mer. Measured by brute force over
# the 29,903-base reference, scoring every oligo at every offset that is not its own
# locus, the best WRONG site is:
#
#     N2-R   (18)  5 mismatches  28%        N3-F   (22)  6 mismatches  27%
#     N2-F   (20)  6 mismatches  30%        RdRp-P (25)  8 mismatches  32%
#     E-P    (26) 11 mismatches  42%        N1-P   (24) 10 mismatches  42%
#
# Coincidence starts at 27% of the oligo's length, and the flat cap of 8 admitted a
# coincidental site for 10 of the 15 oligos. That is not hypothetical: masking all 18
# bases of N2-R's site with N made `find` report position 20,776 - 8,436 bases from the
# real site - with 7 mismatches and `usable` true, instead of saying the site could not be
# read. The ambiguous alignment at the true locus was already being kept as a fallback and
# lost the ranking to the coincidence, because `_preference` adds unknowns to mismatches
# and 7 beats 18.
#
# A fifth of the oligo leaves clear air under every measured coincidence while staying
# well above real substitution drift, which is a handful of bases at most - a primer
# diverging by a fifth would not amplify anyway.
MAX_MISMATCH_SHARE = 0.20
MIN_MISMATCH_CAP = 2


def mismatch_cap(oligo_length: int) -> int:
    """How many definite mismatches a site may have and still be this oligo's site."""
    return min(MAX_MISMATCHES, max(MIN_MISMATCH_CAP, int(oligo_length * MAX_MISMATCH_SHARE)))


# Indels an alignment may use to explain a site before giving up. Three covers the
# in-frame deletions that actually circulate (Alpha's ORF1a delta-3675-3677, S
# delta-69-70); beyond that the site is not this primer's binding site any more.
MAX_GAPS = 3

# Affine gap costs, against a mismatch's 1: opening a gap costs GAP_OPEN and each
# further base GAP_EXTEND, so one three-base deletion is one event (2.0 + 2 x 0.5 = 3.0)
# rather than three. A flat per-base cost charged it 6.0, exactly tying the six
# mismatches it explained, and the alignment lost on the tie.
GAP_OPEN = 2.0
GAP_EXTEND = 0.5


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Hit:
    """The best place this oligo binds, and how well."""

    position: int
    mismatches: int  # definite: target base is known and disagrees
    ambiguous: int  # target base is N - unknown, counted, never guessed
    three_prime_mismatches: int
    found: bool = True
    role: str = ""
    # Indels the alignment needed to explain this site. Zero on the ordinary
    # substitution path; non-zero means the site was explained by a deletion or
    # insertion rather than by a run of mismatches.
    indels: int = 0

    @property
    def usable(self) -> bool:
        """Whether this measurement can go into a rate.

        An alignment with unknown bases under it is not evidence either way.
        """
        return self.found and self.ambiguous == 0

    @property
    def blocks_extension(self) -> int:
        """3' mismatches that can actually stop the reaction.

        **Only on a primer.** The 3'-end rule exists because polymerase extends
        from that terminus, so a mismatch there prevents extension. A hydrolysis
        probe is never extended - its 3' end carries the quencher and is
        chemically blocked - so a mismatch there weakens binding like any other
        and does not stop amplification.

        Counting probe 3' mismatches as assay-breaking is what first made the
        CDC N1 assay look like it had stopped working outright in 2022: the
        Omicron mutation under its probe sits two bases from the probe's 3' end,
        which is meaningless for a probe and fatal for a primer.
        """
        return 0 if self.role == "probe" else self.three_prime_mismatches


NOT_FOUND = Hit(
    position=-1,
    mismatches=MAX_MISMATCHES + 1,
    ambiguous=0,
    three_prime_mismatches=THREE_PRIME_WINDOW,
    found=False,
)


def _pattern(oligo: str) -> re.Pattern[str]:
    """Regex matching this oligo exactly, treating target `N` as a wildcard.

    `N` is added to every character class so a sequencing gap does not prevent
    the binding site from being *located*. Whether those positions matched is a
    separate question, answered by `score_at` and kept out of the rates.
    """
    return re.compile("".join(f"[{IUPAC[base]}N]" for base in oligo.upper()))


# Mismatches the seed search is guaranteed to see through. With m mismatches an oligo
# splits into m+1 clean runs totalling len-m bases, so the longest is at least
# (len - m) // (m + 1) - and a seed no longer than that must land inside one of them.
SEED_GUARANTEE = 2
MIN_SEED = 5
MAX_SEED = 8


def _seed_length(oligo_length: int) -> int:
    """Longest seed that still cannot be missed by SEED_GUARANTEE mismatches."""
    safe = (oligo_length - SEED_GUARANTEE) // (SEED_GUARANTEE + 1)
    return max(MIN_SEED, min(MAX_SEED, safe))


def _seeds(oligo: str, k: int | None = None) -> list[tuple[int, str]]:
    """Unambiguous k-mers of the oligo, with their offsets.

    Used to find candidate positions fast: searching seeds with `str.find`, which runs at
    C speed, beats sliding a window over 30,000 bases in Python by a wide margin.

    The length and the stride both matter, and both used to be wrong. Fixed 8-mers at a
    stride of 4 gave an 18-mer exactly three seeds - 0-7, 4-11, 8-15 - which do not tile
    it and leave bases 16-17 inside no seed at all, so two substitutions could break every
    seed and the site came back NOT_FOUND. Enumerated over every pair of positions, that
    was **31.4% of all two-mismatch sites in an 18-mer**, 6.9% in a 22-mer, and 60.8% of
    three-mismatch sites in an 18-mer. CDC N2's reverse primer is an 18-mer, and
    NOT_FOUND is read as "this assay has no binding site" - a far stronger claim than "it
    has two substitutions", and the one this tool exists to tell apart.

    So the seed is sized to the oligo rather than fixed, and every offset is used rather
    than every other one. The docstring's appeal to the pigeonhole principle only holds
    once both are true.
    """
    if k is None:
        k = _seed_length(len(oligo))
    out = []
    for start in range(0, len(oligo) - k + 1):
        piece = oligo[start : start + k]
        if all(base in "ACGT" for base in piece):
            out.append((start, piece))
    return out


def score_at(
    oligo: str, target: str, position: int, three_prime_at_start: bool = False
) -> Hit | None:
    """Score this oligo against the target starting at `position`.

    `three_prime_at_start` says which end of the *search string* is the oligo's
    3' terminus. It is False for an oligo searched as written, and True for one
    searched as its reverse complement: reverse-complementing turns a sequence
    end for end, so the primer's 3' base becomes index 0 of the string being
    searched. Getting this wrong evaluates the extension-blocking rule - the
    whole point of the 3' window - at the harmless end of every reverse primer.
    """
    window = target[position : position + len(oligo)]
    if len(window) != len(oligo):
        return None

    mismatches = ambiguous = three_prime = 0
    last = len(oligo) - 1
    for index, (base, actual) in enumerate(zip(oligo.upper(), window.upper(), strict=True)):
        if actual == "N":
            ambiguous += 1
            continue
        if actual in IUPAC[base]:
            continue
        mismatches += 1
        distance_from_3_prime = index if three_prime_at_start else last - index
        if distance_from_3_prime < THREE_PRIME_WINDOW:
            three_prime += 1
    return Hit(
        position=position,
        mismatches=mismatches,
        ambiguous=ambiguous,
        three_prime_mismatches=three_prime,
    )


def _align_gapped(
    oligo: str, window: str, three_prime_at_start: bool, max_gaps: int = MAX_GAPS
) -> Hit | None:
    """Score this oligo against a window allowing up to `max_gaps` indels.

    `score_at` compares base i of the oligo with base i of the target, which is
    the right thing for a substitution and the wrong thing for an indel: every
    position after a deletion is compared against its neighbour, so one deleted
    base reads as a run of mismatches. A single base deleted in the middle of a
    21-mer was reported as **eight** mismatches with **three** of them inside
    the 3' window - and a 3' mismatch is the one signal this tool treats as
    stopping the reaction. So a real deletion in the middle of a primer site
    was reported as the primer's business end being destroyed.

    Deletions are not hypothetical in this domain: Alpha carried ORF1a
    delta-3675-3677 and S delta-69-70, and a deletion inside a primer site is
    exactly the kind of event an assay owner needs told apart from drift.

    Gap costs are affine: GAP_OPEN to start a gap and GAP_EXTEND per extra base, so
    one three-base deletion is charged as a single event rather than as three. A flat
    per-base cost made a three-base deletion cost exactly as much as the six mismatches
    it explained, so the alignment was rejected on a tie and the deletion stayed
    misreported. A mismatch costs 1, so no gap is opened to save a single substitution
    and substitution-only drift - the common case - stays on the ungapped path.

    Returns None when no gapped alignment consumes the whole oligo, or when it needed
    more than `max_gaps` gap bases.
    """
    n, m = len(oligo), len(window)
    if n == 0 or m == 0:
        return None
    oligo, window = oligo.upper(), window.upper()
    big = float("inf")

    # Gotoh's three states: M aligns a base to a base, D puts the oligo base over a gap
    # (deleted from the target), I puts a target base over a gap (inserted in it).
    # Separate states are what make the gap penalty affine - a flat matrix cannot tell
    # "opening a gap" from "continuing one".
    M = [[big] * (m + 1) for _ in range(n + 1)]
    D = [[big] * (m + 1) for _ in range(n + 1)]
    I = [[big] * (m + 1) for _ in range(n + 1)]  # noqa: E741 - the standard name
    M[0][0] = 0.0
    for j in range(1, m + 1):
        I[0][j] = GAP_OPEN + GAP_EXTEND * (j - 1)
    for i in range(1, n + 1):
        D[i][0] = GAP_OPEN + GAP_EXTEND * (i - 1)

    for i in range(1, n + 1):
        base = oligo[i - 1]
        for j in range(1, m + 1):
            actual = window[j - 1]
            step = 0.0 if (actual == "N" or actual in IUPAC[base]) else 1.0
            M[i][j] = min(M[i - 1][j - 1], D[i - 1][j - 1], I[i - 1][j - 1]) + step
            D[i][j] = min(D[i - 1][j] + GAP_EXTEND, M[i - 1][j] + GAP_OPEN,
                          I[i - 1][j] + GAP_OPEN)
            I[i][j] = min(I[i][j - 1] + GAP_EXTEND, M[i][j - 1] + GAP_OPEN,
                          D[i][j - 1] + GAP_OPEN)

    ends = [(min(M[n][j], D[n][j], I[n][j]), abs(j - n), j)
            for j in range(max(0, n - max_gaps), min(m, n + max_gaps) + 1)]
    ends = [e for e in ends if e[0] < big]
    if not ends:
        return None
    _, _, end_j = min(ends)

    # Traceback through the state that actually produced each cell.
    mismatches = ambiguous = three_prime = gaps = 0
    i, j = n, end_j
    state = min(("M", M[i][j]), ("D", D[i][j]), ("I", I[i][j]), key=lambda kv: kv[1])[0]
    while i > 0 or j > 0:
        if state == "M":
            if i == 0 or j == 0:
                state = "D" if j == 0 else "I"
                continue
            base, actual = oligo[i - 1], window[j - 1]
            if actual == "N":
                ambiguous += 1
            elif actual not in IUPAC[base]:
                mismatches += 1
                # Measured on the OLIGO index, which the alignment preserves - that is
                # the whole reason for aligning before counting 3' mismatches.
                distance = (i - 1) if three_prime_at_start else (n - i)
                if distance < THREE_PRIME_WINDOW:
                    three_prime += 1
            step = 0.0 if (actual == "N" or actual in IUPAC[base]) else 1.0
            previous = M[i][j] - step
            i, j = i - 1, j - 1
            state = min(("M", M[i][j]), ("D", D[i][j]), ("I", I[i][j]),
                        key=lambda kv: abs(kv[1] - previous))[0]
        elif state == "D":
            gaps += 1
            came_from_extend = i > 1 and abs(D[i][j] - (D[i - 1][j] + GAP_EXTEND)) < 1e-9
            i -= 1
            state = "D" if came_from_extend else "M"
            if i == 0 and j == 0:
                break
        else:
            gaps += 1
            came_from_extend = j > 1 and abs(I[i][j] - (I[i][j - 1] + GAP_EXTEND)) < 1e-9
            j -= 1
            state = "I" if came_from_extend else "M"
            if i == 0 and j == 0:
                break
    if gaps == 0 or gaps > max_gaps:
        return None
    return Hit(
        position=-1,
        mismatches=mismatches,
        ambiguous=ambiguous,
        three_prime_mismatches=three_prime,
        indels=gaps,
    )


def _preference(hit: Hit) -> tuple[int, int]:
    """How good an alignment is, lower being better.

    Fewest *unconfirmed* positions first, then fewest real mismatches. A
    position is confirmed only when the target says a definite base and that
    base is one the oligo can be; both a mismatch and an `N` leave it
    unconfirmed.

    Ranking on mismatches alone is what let a run of `N` beat the real binding
    site: `_pattern` treats target `N` as a wildcard, so a gap of `N` as long as
    the oligo scores zero mismatches anywhere in the genome. Counting unknowns
    against an alignment is what distinguishes "this site matches" from "this
    site says nothing".
    """
    return (hit.mismatches + hit.ambiguous, hit.mismatches)


def find(oligo: str, target: str, three_prime_at_start: bool = False) -> Hit:
    """Best binding site for this oligo on this strand of the target.

    `three_prime_at_start` is passed through to `score_at`; see there. Callers
    that have already reverse-complemented a reverse primer must set it, and
    `find_oligo` does it for them.
    """
    oligo = oligo.upper()
    target = target.upper()
    bad = set(oligo) - set(IUPAC)
    if bad:
        raise ValueError(f"not IUPAC nucleotide codes: {sorted(bad)}")
    if not oligo or len(target) < len(oligo):
        return NOT_FOUND

    best: Hit | None = None
    cap = mismatch_cap(len(oligo))

    # Fast path: a perfect site, allowing N. Most genomes, most assays.
    exact = _pattern(oligo).search(target)
    if exact is not None:
        scored = score_at(oligo, target, exact.start(), three_prime_at_start)
        if scored is not None:
            if scored.mismatches == 0 and scored.ambiguous == 0:
                return scored
            # Zero mismatches over unknown bases is not a match, it is no
            # evidence - and `re.search` returns the leftmost hit, so a
            # sequencing gap early in the genome would otherwise be preferred
            # over the real binding site. Keep it only as a fallback for a
            # genome that has nothing better anywhere.
            best = scored

    seen: set[int] = set()
    # Positions a seed pointed at that the mismatch cap rejected. An indel inside the
    # site puts the RIGHT position in here, so they are kept for the gapped pass.
    rejected: set[int] = set()
    for offset, seed in _seeds(oligo):
        start = target.find(seed)
        while start != -1:
            position = start - offset
            if position >= 0 and position not in seen:
                seen.add(position)
                scored = score_at(oligo, target, position, three_prime_at_start)
                # Anything above the cap is noise found elsewhere in a 30,000
                # base genome, so it is dropped here rather than after the
                # ranking: a coincidence has few unknowns and would otherwise
                # outrank a real but gappy or unreadable site, turning "this
                # genome cannot say" into "this assay has no binding site".
                # The cap is a share of the oligo's length, not a flat 8 - see
                # MAX_MISMATCH_SHARE for the measurement that fixed it there.
                if scored is None or scored.mismatches > cap:
                    # Keep the position anyway. An indel inside the site makes every
                    # base after it compare against its neighbour, so the ungapped
                    # score at the RIGHT place is terrible and gets dropped here - and
                    # then nothing is left and the site reports as absent. A 3-base
                    # deletion in a 22-mer primer did exactly that: NOT_FOUND, which
                    # reads as "this assay has no binding site" rather than "there is a
                    # deletion in it". The gapped pass below is the only thing that can
                    # tell those apart, and it has to be allowed to see these.
                    if position >= 0:
                        rejected.add(position)
                    start = target.find(seed, start + 1)
                    continue
                if best is None or _preference(scored) < _preference(best):
                    best = scored
            start = target.find(seed, start + 1)

    if best is None and not rejected:
        return NOT_FOUND

    # Only now, and only if the best ungapped alignment looks bad enough that an indel
    # is the likelier explanation. A site with two or fewer mismatches is ordinary
    # substitution drift and must not be refitted as a gap; above that, one deleted base
    # can masquerade as a run of mismatches and - worse - manufacture mismatches inside
    # the 3' window, which is the signal this tool reports as stopping the reaction.
    # Rejected positions are a last resort, used only when nothing survived the cap.
    # Letting them compete with a surviving alignment lets a gapped NOISE site beat a
    # genuinely unreadable one: a binding site that is entirely `N` is "cannot say", and
    # promoting a 5-mismatch-plus-a-gap coincidence over it turns that into "the assay
    # has lost its target", which is the opposite conclusion.
    anchors: set[int] = set()
    if best is None:
        anchors |= rejected
    elif best.position >= 0 and best.mismatches > 2:
        anchors.add(best.position)
    if anchors:
        # The start has to be searched, not assumed. An indel shifts everything after it,
        # so the best UNGAPPED alignment of a site containing one lands off the true
        # start - a 3-base deletion put it 3 bases early and a 3-base insertion 3 late -
        # and a window anchored on that start cannot be repaired by gaps alone.
        candidates = []
        for anchor in anchors:
            for shift in range(-MAX_GAPS, MAX_GAPS + 1):
                start = anchor + shift
                if start < 0:
                    continue
                window = target[start : start + len(oligo) + MAX_GAPS]
                if len(window) < len(oligo) - MAX_GAPS:
                    continue
                gapped = _align_gapped(oligo, window, three_prime_at_start)
                if gapped is None:
                    continue
                # Affine, matching _align_gapped: one run of gaps is one event. Scoring
                # per gap base here would reject a 3-base deletion that the aligner had
                # already decided was the better explanation.
                score = (
                    gapped.mismatches
                    + gapped.ambiguous
                    + GAP_OPEN
                    + GAP_EXTEND * max(0, gapped.indels - 1)
                )
                candidates.append((score, abs(shift), start, gapped))
        if candidates:
            score, _, start, gapped = min(candidates, key=lambda c: (c[0], c[1]))
            # Against the ungapped alternative, or against the cap when there was no
            # usable ungapped alignment at all.
            ceiling = (
                best.mismatches + best.ambiguous if best is not None else cap + 1
            )
            if score < ceiling:
                return replace(gapped, position=start)
    if best is None:
        return NOT_FOUND
    return best


def find_oligo(sequence: str, oligo: str, role: str) -> Hit:
    """Locate an oligo, taking the strand from its role.

    A reverse primer binds the template strand, so it is its reverse complement
    that appears in a forward-strand genome record. Searching for it as written
    finds nothing and reports every assay as failed.

    Reverse-complementing also flips which end of the search string is the
    primer's 3' terminus, which is why `three_prime_at_start` travels with it:
    without that, the extension-blocking rule is applied to the 5' end of every
    reverse primer, so a real 3'-terminal mismatch scores zero and a harmless
    5' one is reported as assay-breaking.
    """
    reverse = role == "reverse"
    search_for = reverse_complement(oligo) if reverse else oligo
    hit = find(search_for, sequence, three_prime_at_start=reverse)
    # The role travels with the hit, because whether a 3' mismatch matters
    # depends on it: a primer gets extended from that end, a probe does not.
    return replace(hit, role=role)
