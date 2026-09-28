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
MAX_MISMATCHES = 8


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


def _seeds(oligo: str, k: int = 8) -> list[tuple[int, str]]:
    """Unambiguous k-mers of the oligo, with their offsets.

    Used to find candidate positions fast. By the pigeonhole principle an
    alignment with few mismatches must contain at least one exact seed, so
    searching seeds with `str.find` (which runs at C speed) beats sliding a
    window over 30,000 bases in Python by a wide margin.
    """
    out = []
    for start in range(0, len(oligo) - k + 1, max(1, k // 2)):
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
    for offset, seed in _seeds(oligo):
        start = target.find(seed)
        while start != -1:
            position = start - offset
            if position >= 0 and position not in seen:
                seen.add(position)
                scored = score_at(oligo, target, position, three_prime_at_start)
                # Anything above MAX_MISMATCHES is noise found elsewhere in a
                # 30,000 base genome, so it is dropped here rather than after
                # the ranking: a 9-mismatch coincidence has few unknowns and
                # would otherwise outrank a real but gappy site, turning "this
                # genome cannot say" into "this assay has no binding site".
                if scored is None or scored.mismatches > MAX_MISMATCHES:
                    start = target.find(seed, start + 1)
                    continue
                if best is None or _preference(scored) < _preference(best):
                    best = scored
            start = target.find(seed, start + 1)

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
