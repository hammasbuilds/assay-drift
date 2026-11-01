"""Nearest-neighbour duplex thermodynamics, for primers against what is actually there.

Section 4 of the paper is about `likely_failing` resting on a positional rule - a mismatch
inside the last five bases of a primer's 3' end counts as blocking - whose threshold is a
primer-design rule of thumb rather than a calibrated quantity. The honest response to "that
number is uncalibrated" is not to pick a better number by judgement. It is to compute
something that is not a judgement at all and report where the two disagree.

The nearest-neighbour model is that quantity. A duplex's stability is the sum of
contributions from each overlapping pair of stacked base pairs, plus initiation terms, and
the parameters are measured by UV melting rather than chosen:

* Watson-Crick stacks - Allawi & SantaLucia (1997), *Biochemistry* 36:10581.
* Internal mismatches - Allawi & SantaLucia (1997-1998), Peyret et al. (1999),
  Watkins & SantaLucia (2005).
* Terminal mismatches - SantaLucia & Peyret (2001).
* Dangling ends - Bommarito et al. (2000).
* Salt correction and the unified framework - SantaLucia (1998), *PNAS* 95:1460, and
  SantaLucia & Hicks (2004), *Annu Rev Biophys Biomol Struct* 33:415.

The quantity this module exists to produce is `destabilisation`: the free-energy cost, in
kcal/mol, of the observed target sequence relative to a perfectly complementary one.
dGG is a difference between two numbers from the same model, so the parameter set's own
error - the mismatch papers report about 6.4% average deviation on dG37 - largely cancels,
which is why the difference is reported and the absolute dG is not.

What this is not: a prediction of whether a PCR reaction fails. Amplification depends on
polymerase, cycling conditions and extension from the 3' terminus, none of which is in a
two-state duplex model. dGG ranks mismatches by how much binding they cost. The paper uses
it as a second, independent opinion on the positional rule, and says so.

The tables are vendored rather than imported, because this package has no runtime
dependencies. They were dumped from Biopython's copies rather than retyped from the
papers - a wrong digit in a hand-transcribed table is a defect nothing would ever find -
and `tests/test_thermo.py` computes every quantity here against Biopython to prove the
walk over them agrees with an independent implementation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

# The gas constant in cal/(K*mol), as the source papers use it.
R: Final = 1.987

# Standard-state reference temperature for dG, in degrees Celsius. Every dG and dGG this
# module reports is dG at this temperature unless the caller says otherwise. 37 degrees is
# the convention of the source papers; a PCR annealing step is nearer 55-60, and
# `destabilisation` takes a temperature for that reason.
REFERENCE_TEMPERATURE_C: Final = 37.0

# The reaction the published assays specify. Both CDC and Charite protocols are ordinary
# TaqMan chemistry; these are the concentrations their inserts give, in millimolar for the
# ions and nanomolar for the oligos.
DEFAULT_NA_MM: Final = 50.0
DEFAULT_MG_MM: Final = 3.0
DEFAULT_DNTP_MM: Final = 0.8
DEFAULT_OLIGO_NM: Final = 400.0

# Allawi & SantaLucia (1997), Biochemistry 36:10581 - Watson-Crick
# Biopython 1.88 `MeltingTemp.DNA_NN3`, 17 entries.
NN: Final[dict[str, tuple[float, float]]] = {
    "AA/TT": (-7.9, -22.2),
    "AT/TA": (-7.2, -20.4),
    "CA/GT": (-8.5, -22.7),
    "CG/GC": (-10.6, -27.2),
    "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2),
    "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9),
    "GT/CA": (-8.4, -22.4),
    "TA/AT": (-7.2, -21.3),
    "init": (0, 0),
    "init_5T/A": (0, 0),
    "init_A/T": (2.3, 4.1),
    "init_G/C": (0.1, -2.8),
    "init_allA/T": (0, 0),
    "init_oneG/C": (0, 0),
    "sym": (0, -1.4),
}

# Allawi & SantaLucia (1997-1998); Peyret et al. (1999); Watkins & SantaLucia (2005)
# Biopython 1.88 `MeltingTemp.DNA_IMM1`, 87 entries.
INTERNAL_MISMATCH: Final[dict[str, tuple[float, float]]] = {
    "AA/TA": (1.2, 1.7),
    "AA/TC": (2.3, 4.6),
    "AA/TG": (-0.6, -2.3),
    "AA/TI": (-0.7, -2.6),
    "AC/TA": (5.3, 14.6),
    "AC/TC": (0.0, -4.4),
    "AC/TI": (-8.8, -25.4),
    "AC/TT": (0.7, 0.2),
    "AG/TA": (-0.7, -2.3),
    "AG/TG": (-3.1, -9.5),
    "AG/TI": (0.1, -1.8),
    "AG/TT": (1.0, 0.9),
    "AI/TA": (-8.3, -25.0),
    "AI/TC": (-8.9, -25.5),
    "AI/TG": (-4.9, -15.8),
    "AI/TI": (-3.3, -11.9),
    "AI/TT": (0.49, -0.7),
    "AT/TC": (-1.2, -6.2),
    "AT/TG": (-2.5, -8.3),
    "AT/TI": (-5.6, -18.7),
    "AT/TT": (-2.7, -10.8),
    "CA/GA": (-0.9, -4.2),
    "CA/GC": (1.9, 3.7),
    "CA/GG": (-0.7, -2.3),
    "CA/GI": (-7.0, -20.0),
    "CC/GA": (0.6, -0.6),
    "CC/GC": (-1.5, -7.2),
    "CC/GI": (-8.3, -23.8),
    "CC/GT": (-0.8, -4.5),
    "CG/GA": (-4.0, -13.2),
    "CG/GG": (-4.9, -15.3),
    "CG/GI": (5.8, 16.9),
    "CG/GT": (-4.1, -11.7),
    "CI/GA": (2.6, 8.9),
    "CI/GC": (-5.4, -13.7),
    "CI/GG": (7.1, 21.3),
    "CI/GI": (1.3, 3.0),
    "CI/GT": (-1.0, -2.4),
    "CT/GC": (-1.5, -6.1),
    "CT/GG": (-2.8, -8.0),
    "CT/GI": (0.1, -1.0),
    "CT/GT": (-5.0, -15.8),
    "GA/CA": (-2.9, -9.8),
    "GA/CC": (5.2, 14.2),
    "GA/CG": (-0.6, -1.0),
    "GA/CI": (-7.6, -20.2),
    "GC/CA": (-0.7, -3.8),
    "GC/CC": (3.6, 8.9),
    "GC/CI": (-5.0, -12.6),
    "GC/CT": (2.3, 5.4),
    "GG/CA": (0.5, 3.2),
    "GG/CG": (-6.0, -15.8),
    "GG/CI": (-7.6, -22.0),
    "GG/CT": (3.3, 10.4),
    "GG/TT": (5.8, 16.3),
    "GI/CA": (-7.8, -21.1),
    "GI/CC": (-6.8, -19.1),
    "GI/CG": (-1.1, -3.2),
    "GI/CI": (-0.5, -1.3),
    "GI/CT": (-3.5, -10.6),
    "GT/CC": (5.2, 13.5),
    "GT/CG": (-4.4, -12.3),
    "GT/CI": (-4.3, -12.1),
    "GT/CT": (-2.2, -8.4),
    "GT/TG": (4.1, 9.5),
    "TA/AA": (4.7, 12.9),
    "TA/AC": (3.4, 8.0),
    "TA/AG": (0.7, 0.7),
    "TA/AI": (-1.3, -4.6),
    "TC/AA": (7.6, 20.2),
    "TC/AC": (6.1, 16.4),
    "TC/AI": (-4.9, -13.9),
    "TC/AT": (1.2, 0.7),
    "TG/AA": (3.0, 7.4),
    "TG/AG": (1.6, 3.6),
    "TG/AI": (1.0, 1.0),
    "TG/AT": (-0.1, -1.7),
    "TG/GT": (-1.4, -6.2),
    "TI/AA": (-3.4, -11.2),
    "TI/AC": (-5.9, -17.4),
    "TI/AG": (-1.9, -8.5),
    "TI/AI": (0.1, -2.3),
    "TI/AT": (-6.5, -22.0),
    "TT/AC": (1.0, 0.7),
    "TT/AG": (-1.3, -5.3),
    "TT/AI": (-0.8, -4.3),
    "TT/AT": (0.2, -1.5),
}

# SantaLucia & Peyret (2001), patent WO 01/94611
# Biopython 1.88 `MeltingTemp.DNA_TMM1`, 48 entries.
TERMINAL_MISMATCH: Final[dict[str, tuple[float, float]]] = {
    "AA/TA": (-3.1, -7.8),
    "AA/TC": (-1.6, -4.0),
    "AA/TG": (-1.9, -4.4),
    "AC/TA": (-1.8, -3.8),
    "AC/TC": (-0.1, 0.5),
    "AC/TT": (-0.9, -1.7),
    "AG/TA": (-2.5, -5.9),
    "AG/TG": (-1.1, -2.1),
    "AG/TT": (-3.2, -8.7),
    "AT/TC": (-2.3, -6.3),
    "AT/TG": (-3.5, -9.4),
    "AT/TT": (-2.4, -6.5),
    "CA/GA": (-4.3, -10.7),
    "CA/GC": (-2.6, -5.9),
    "CA/GG": (-3.9, -9.6),
    "CC/GA": (-2.7, -6.0),
    "CC/GC": (-2.1, -5.1),
    "CC/GT": (-3.2, -8.0),
    "CG/GA": (-6.0, -15.5),
    "CG/GG": (-3.8, -9.5),
    "CG/GT": (-3.8, -9.0),
    "CT/GC": (-3.9, -10.6),
    "CT/GG": (-6.6, -18.7),
    "CT/GT": (-6.1, -16.9),
    "GA/CA": (-8.0, -22.5),
    "GA/CC": (-5.0, -13.8),
    "GA/CG": (-4.3, -11.1),
    "GC/CA": (-3.2, -7.1),
    "GC/CC": (-3.9, -10.6),
    "GC/CT": (-4.9, -13.5),
    "GG/CA": (-4.6, -11.4),
    "GG/CG": (-0.7, -19.2),
    "GG/CT": (-5.7, -15.9),
    "GT/CC": (-3.0, -7.8),
    "GT/CG": (-5.9, -16.1),
    "GT/CT": (-7.4, -21.2),
    "TA/AA": (-2.5, -6.3),
    "TA/AC": (-2.3, -5.9),
    "TA/AG": (-2.0, -4.7),
    "TC/AA": (-2.7, -7.0),
    "TC/AC": (-0.7, -1.3),
    "TC/AT": (-2.5, -6.3),
    "TG/AA": (-2.4, -5.8),
    "TG/AG": (-1.1, -2.7),
    "TG/AT": (-3.9, -10.5),
    "TT/AC": (-0.7, -1.2),
    "TT/AG": (-3.6, -9.8),
    "TT/AT": (-3.2, -8.9),
}


# Degenerate positions in a published primer are a mixture the manufacturer synthesises,
# not an unknown: the oligo pool contains every resolution. Enumerating them is only safe
# while the pool is small, and the catalogue's worst oligo has three degenerate positions.
MAX_RESOLUTIONS: Final = 64

_COMPLEMENT: Final = {"A": "T", "T": "A", "G": "C", "C": "G"}

# Exactly the four bases the parameter tables cover. `N`, and every other IUPAC code, is
# a refusal rather than a guess: see `duplex`.
_CONCRETE: Final = frozenset("ACGT")

# The degenerate codes a published primer may legitimately contain. Deliberately not the
# full IUPAC table from `match`: that one also has to read whatever a GenBank record
# contains, and this one only ever reads a catalogue entry.
_DEGENERATE: Final[dict[str, str]] = {
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
}


def complement(sequence: str) -> str | None:
    """The complementary bases, column for column, with no reversal.

    The tables are keyed on a duplex written as top strand 5' to 3' over bottom strand
    3' to 5', which is the same column order. A reversal here would silently score every
    duplex backwards, which for a palindromic stack would even give the right answer
    some of the time.
    """
    out = []
    for base in sequence.upper():
        mapped = _COMPLEMENT.get(base)
        if mapped is None:
            return None
        out.append(mapped)
    return "".join(out)


def resolutions(oligo: str) -> list[str] | None:
    """Every concrete sequence the synthesised pool for `oligo` contains.

    `None` when the oligo holds a code that is not a nucleotide at all, or when the pool
    is larger than `MAX_RESOLUTIONS`. Returning `None` rather than a truncated list
    matters: a truncated pool would quietly report the best of an arbitrary subset.
    """
    pools = []
    for base in oligo.upper():
        if base in _CONCRETE:
            pools.append(base)
        elif base in _DEGENERATE:
            pools.append(_DEGENERATE[base])
        else:
            return None
    size = 1
    for pool in pools:
        size *= len(pool)
        if size > MAX_RESOLUTIONS:
            return None
    out = [""]
    for pool in pools:
        out = [prefix + base for prefix in out for base in pool]
    return out


def duplex(oligo: str, template: str) -> tuple[float, float] | None:
    """dH in kcal/mol and dS in cal/(K*mol) for one concrete duplex.

    `oligo` is the primer or probe 5' to 3'. `template` is the strand it anneals to,
    written 3' to 5' so that column *i* of each string is one base pair - which is what
    `complement` produces from a same-sense genome slice.

    `None` when either strand holds anything but A, C, G or T. The caller's job is to
    have resolved the oligo's degeneracies and to have refused a target with unknown
    bases; an `N` scored as some particular base is a fabricated measurement, and this
    model has no parameter for one.

    `None` too for a *tandem* mismatch - two adjacent mispaired positions. The published
    parameters cover a mismatch flanked by matched pairs and do not cover two in a row,
    so such a duplex is not measurable with this model rather than measurable at some
    cost. See the refusal at the end of the zipping loop.

    The walk is the published one: strip and score a terminal mismatch at each end, add
    the initiation terms, then sum the stacks. Where a stack is a mismatch the internal
    mismatch table supplies it, in either orientation, since a duplex read from the other
    end is the same duplex.
    """
    top, bottom = oligo.upper(), template.upper()
    if len(top) != len(bottom) or len(top) < 2:
        return None
    if not (set(top) <= _CONCRETE and set(bottom) <= _CONCRETE):
        return None

    delta_h = 0.0
    delta_s = 0.0

    # A mismatched pair at either terminus is not a stack like the others - it frays -
    # and has its own parameter set. Scored, then removed before the zipping.
    left = bottom[:2][::-1] + "/" + top[:2][::-1]
    if left in TERMINAL_MISMATCH:
        delta_h += TERMINAL_MISMATCH[left][0]
        delta_s += TERMINAL_MISMATCH[left][1]
        top, bottom = top[1:], bottom[1:]
    right = top[-2:] + "/" + bottom[-2:]
    if right in TERMINAL_MISMATCH:
        delta_h += TERMINAL_MISMATCH[right][0]
        delta_s += TERMINAL_MISMATCH[right][1]
        top, bottom = top[:-1], bottom[:-1]

    # Initiation. Read off the whole oligo, not the trimmed duplex, as the papers define
    # it and as the reference implementation computes it.
    whole = oligo.upper()
    delta_h += NN["init"][0]
    delta_s += NN["init"][1]
    key = "init_allA/T" if set(whole) <= {"A", "T"} else "init_oneG/C"
    delta_h += NN[key][0]
    delta_s += NN[key][1]
    if whole.startswith("T"):
        delta_h += NN["init_5T/A"][0]
        delta_s += NN["init_5T/A"][1]
    if whole.endswith("A"):
        delta_h += NN["init_5T/A"][0]
        delta_s += NN["init_5T/A"][1]
    ends = whole[0] + whole[-1]
    at = ends.count("A") + ends.count("T")
    gc = ends.count("G") + ends.count("C")
    delta_h += NN["init_A/T"][0] * at + NN["init_G/C"][0] * gc
    delta_s += NN["init_A/T"][1] * at + NN["init_G/C"][1] * gc

    for index in range(len(top) - 1):
        stack = top[index : index + 2] + "/" + bottom[index : index + 2]
        for table in (INTERNAL_MISMATCH, NN):
            if stack in table:
                delta_h += table[stack][0]
                delta_s += table[stack][1]
                break
            if stack[::-1] in table:
                delta_h += table[stack[::-1]][0]
                delta_s += table[stack[::-1]][1]
                break
        else:
            # Reached, and not rare: the published parameter set has no values for a
            # tandem mismatch - two adjacent positions both mispaired, which makes a
            # stack like `TA/CA` that no table covers. Nearest-neighbour additivity is
            # not claimed to hold there, and the reference implementation raises on the
            # same input rather than inventing a value.
            #
            # This is a refusal, not a zero. A silent zero would read as a free
            # mismatch, which is the worst available direction for the error: the
            # duplex the model cannot score at all would come out looking undamaged.
            # Callers must treat `None` as "not measurable", and the corpus does
            # contain the case - C28311T and C28312T are adjacent bases of the same
            # N1 probe.
            return None

    return delta_h, delta_s


def salt_correction(oligo: str, na_mm: float, mg_mm: float, dntp_mm: float) -> float:
    """The entropy correction for the reaction's ions, in cal/(K*mol).

    SantaLucia (1998) equation 5, the length-dependent form, with divalent magnesium
    folded into a sodium equivalent by von Ahsen et al. (2001). dNTPs chelate magnesium,
    so only the excess counts; at or below the dNTP concentration the free magnesium is
    treated as irrelevant, which is what the source does.
    """
    monovalent = na_mm
    if mg_mm > dntp_mm:
        monovalent += 120.0 * math.sqrt(mg_mm - dntp_mm)
    if monovalent <= 0:
        raise ValueError("a salt correction needs a non-zero monovalent concentration")
    return 0.368 * (len(oligo) - 1) * math.log(monovalent * 1e-3)


def free_energy(delta_h: float, delta_s: float, celsius: float = REFERENCE_TEMPERATURE_C) -> float:
    """dG at `celsius` in kcal/mol, from dH in kcal/mol and dS in cal/(K*mol)."""
    return delta_h - (celsius + 273.15) * delta_s / 1000.0


def melting_temperature(
    oligo: str,
    template: str,
    *,
    na_mm: float = DEFAULT_NA_MM,
    mg_mm: float = DEFAULT_MG_MM,
    dntp_mm: float = DEFAULT_DNTP_MM,
    oligo_nm: float = DEFAULT_OLIGO_NM,
) -> float | None:
    """Tm in degrees Celsius for one concrete duplex, or `None` if it cannot be scored.

    The two strands are at wildly different concentrations in a PCR - primer in excess,
    template vanishingly rare - so the effective concentration is the primer's, which is
    what the reference implementation's `dnac1 - dnac2/2` reduces to when the template
    term is dropped.
    """
    scored = duplex(oligo, template)
    if scored is None:
        return None
    delta_h, delta_s = scored
    delta_s += salt_correction(oligo, na_mm=na_mm, mg_mm=mg_mm, dntp_mm=dntp_mm)
    return (1000.0 * delta_h) / (delta_s + R * math.log(oligo_nm * 1e-9)) - 273.15


@dataclass(frozen=True)
class Destabilisation:
    """What the observed target costs this oligo, relative to a perfect match.

    `delta_delta_g` is the headline: kcal/mol of binding energy given up, positive when
    the observed sequence binds worse than a perfect complement. `delta_tm` is the same
    fact in degrees, which is the unit a bench scientist reads.

    Both are differences between two numbers out of one parameter set, which is the
    reason to report them rather than either absolute dG. The mismatch parameters carry
    about 6.4% average deviation on dG37; a difference between two duplexes that share
    every stack but one cancels most of that, while an absolute dG does not.

    `resolution` names which member of a degenerate oligo's pool was scored - the one
    that binds the observed target best, because that is the species that amplifies it.
    """

    delta_delta_g: float
    delta_tm: float
    perfect_g: float
    observed_g: float
    perfect_tm: float
    observed_tm: float
    resolution: str

    @property
    def costs_nothing(self) -> bool:
        """Whether the observed target is at least as stable as a perfect complement.

        Not a rounding of `delta_delta_g <= 0`: a degenerate oligo's pool can contain a
        member that matches the observed sequence exactly, in which case the cost really
        is zero and the mismatch the aligner counted is against a different pool member.
        """
        return self.delta_delta_g <= 0.0


def destabilisation(
    oligo: str,
    observed: str,
    *,
    celsius: float = REFERENCE_TEMPERATURE_C,
    na_mm: float = DEFAULT_NA_MM,
    mg_mm: float = DEFAULT_MG_MM,
    dntp_mm: float = DEFAULT_DNTP_MM,
    oligo_nm: float = DEFAULT_OLIGO_NM,
) -> Destabilisation | None:
    """The thermodynamic cost of `observed` to `oligo`, or `None` if it cannot be scored.

    `observed` is the genome slice under the oligo's binding site, in the oligo's own
    sense - which is what this package's aligner reports, since it compares an oligo
    against the strand the oligo's sequence appears in. The strand the oligo actually
    anneals to is that slice's complement, and taking it is this function's job rather
    than the caller's, because getting it backwards is a silent error.

    `None` rather than a number whenever the measurement cannot be made: a target base
    that is not A, C, G or T, an oligo carrying a code that is not a nucleotide, a
    degenerate pool too large to enumerate, or a length disagreement. Every one of those
    is a case where a returned float would be a fabrication, and the corpus contains all
    of them.
    """
    if len(oligo) != len(observed):
        return None
    observed = observed.upper()
    if not set(observed) <= _CONCRETE:
        return None
    pool = resolutions(oligo)
    if not pool:
        return None

    observed_template = complement(observed)
    if observed_template is None:
        return None

    # A pool member identical to what is there has no mismatch against it at all, so the
    # cost is zero by definition and no search is needed. This is not an optimisation: a
    # G-T wobble can be more stable than an A-T pair, so the "best binding" search below
    # can genuinely prefer a *different* species over the one that matches exactly, and
    # would then report that species' own bookkeeping as a cost. Reporting a cost for a
    # target the oligo pool contains an exact complement for would be wrong.
    if observed in pool:
        perfect = duplex(observed, observed_template)
        if perfect is not None:
            observed_tm = melting_temperature(
                observed,
                observed_template,
                na_mm=na_mm,
                mg_mm=mg_mm,
                dntp_mm=dntp_mm,
                oligo_nm=oligo_nm,
            )
            if observed_tm is not None:
                energy = free_energy(*perfect, celsius=celsius)
                return Destabilisation(
                    delta_delta_g=0.0,
                    delta_tm=0.0,
                    perfect_g=energy,
                    observed_g=energy,
                    perfect_tm=observed_tm,
                    observed_tm=observed_tm,
                    resolution=observed,
                )

    best: Destabilisation | None = None
    for candidate in pool:
        observed_pair = duplex(candidate, observed_template)
        perfect_template = complement(candidate)
        if observed_pair is None or perfect_template is None:
            continue
        perfect_pair = duplex(candidate, perfect_template)
        if perfect_pair is None:
            continue
        observed_g = free_energy(*observed_pair, celsius=celsius)
        perfect_g = free_energy(*perfect_pair, celsius=celsius)
        if best is not None and observed_g >= best.observed_g:
            continue
        observed_tm = melting_temperature(
            candidate,
            observed_template,
            na_mm=na_mm,
            mg_mm=mg_mm,
            dntp_mm=dntp_mm,
            oligo_nm=oligo_nm,
        )
        perfect_tm = melting_temperature(
            candidate,
            perfect_template,
            na_mm=na_mm,
            mg_mm=mg_mm,
            dntp_mm=dntp_mm,
            oligo_nm=oligo_nm,
        )
        if observed_tm is None or perfect_tm is None:
            continue
        best = Destabilisation(
            delta_delta_g=observed_g - perfect_g,
            delta_tm=observed_tm - perfect_tm,
            perfect_g=perfect_g,
            observed_g=observed_g,
            perfect_tm=perfect_tm,
            observed_tm=observed_tm,
            resolution=candidate,
        )
    return best
