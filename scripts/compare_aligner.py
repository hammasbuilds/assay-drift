"""Compare this tool's oligo matcher against Biopython's pairwise aligner.

Both answer "where does this oligo sit in this genome, and how many bases disagree". They
answer it differently: `assaydrift.match` searches unambiguous seeds, scores a window with
IUPAC codes read as sets on both sides, and refuses a site whose mismatches exceed a
fifth of the oligo's length; `Bio.Align.PairwiseAligner` runs a full local alignment with
a substitution matrix and affine gaps and returns the best-scoring one wherever it is.

Neither is a baseline for the other, and the point of running them together is not a score.
It is the three things the comparison can actually show:

  * whether they put the oligo in the same place,
  * whether they count the same number of mismatches once they do,
  * and where the purpose-built rules change the answer - an ambiguous subject base, an
    `N` under the oligo, a site too poor to be the binding site at all.

Biopython is not a dependency of this package and never will be. It is installed into its
own virtualenv, `.venv-compare`, which this script creates if it is missing.

    python scripts/compare_aligner.py --genomes 40
    python scripts/compare_aligner.py --genomes 40 --out docs/aligner-comparison.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from assaydrift.match import compare_base, find_oligo, reverse_complement  # noqa: E402
from assaydrift.primers import ASSAYS  # noqa: E402

VENV = ROOT / ".venv-compare"
CACHE = ROOT / "data" / "cache"
SEED = 20261029

# Runs in the other interpreter: Biopython's side of every pair, in one subprocess, so
# the aligner is set up once rather than per oligo.
BIOPYTHON_SIDE = '''
import json, sys, time
from Bio import Align
from Bio.Align import substitution_matrices

blob = json.loads(open(sys.argv[1], encoding="utf-8").read())
# Genomes once, keyed; pairs name one. Sent per pair, 15 oligos meant 15 copies of the
# same 30,000-base record and a payload 15 times the size it needed to be.
targets, work = blob["targets"], blob["pairs"]

# Scoring chosen to ask Biopython the same question, as nearly as its parameters allow:
# a match is +1, a substitution -1, and a gap is expensive enough that the aligner
# prefers a substitution - this tool counts mismatches, and an aligner free to open gaps
# answers a different question. local mode, because an oligo is 20 bases inside 30,000.
def configure(mode):
    a = Align.PairwiseAligner()
    a.mode = mode
    a.match_score = 1
    a.mismatch_score = -1
    # Gaps are made expensive on purpose: this tool counts mismatches, and an aligner
    # free to open gaps answers a different question.
    a.open_gap_score = -10
    a.extend_gap_score = -10
    return a


# Out of the box: local alignment, best-scoring subsequence anywhere in the genome.
local = configure("local")
# Local mode returns the best-scoring SUBSEQUENCE, so a short near-perfect partial match
# outscores a full-length match carrying mismatches: for one primer it reported a site
# 24,800 bases from the real one with 9 of 20 bases unaligned. `unaligned_oligo_bases`
# records that per pair rather than being argued about.
#
# An earlier version ran a second, "anchored" pass meant to force every base to align. It
# returned 19 mismatches of 24 for a probe the matcher scores at 0, because it anchored on
# a position trimming had already shifted and then counted columns of a global alignment
# that is degenerate inside a window of the oligo's own length. The column it was reaching
# for is computed on the other side now, without an aligner: plain base-by-base equality
# at the site this tool reports, which isolates IUPAC awareness and nothing else.

out = []
for item in work:
    oligo, target = item["oligo"], targets[item["accession"]]
    start = time.perf_counter()
    try:
        best = local.align(target, oligo)[0]
    except Exception as exc:
        out.append({"key": item["key"], "error": "%s: %s" % (type(exc).__name__, exc)})
        continue
    elapsed = time.perf_counter() - start
    # The aligned blocks give the target coordinates; mismatches are columns where the
    # letters differ, counted over the aligned region only.
    rows = best.aligned
    mismatches = 0
    aligned_bases = 0
    for (ts, te), (os_, oe) in zip(rows[0], rows[1]):
        for k in range(te - ts):
            aligned_bases += 1
            if target[ts + k] != oligo[os_ + k]:
                mismatches += 1
    position = int(rows[0][0][0]) if len(rows[0]) else None
    out.append({
        "key": item["key"],
        "position": position,
        "mismatches": mismatches,
        "aligned_bases": aligned_bases,
        "unaligned_oligo_bases": len(oligo) - aligned_bases,
        "score": float(best.score),
        "seconds": round(elapsed, 4),
    })

print("__JSON__" + json.dumps(out))
'''


def python_in(venv: Path) -> Path:
    return venv / ("Scripts/python.exe" if os.name == "nt" else "bin/python")


def ensure_biopython() -> str:
    python = python_in(VENV)
    if not python.exists():
        subprocess.run([sys.executable, "-m", "venv", str(VENV)], check=True)
    shown = subprocess.run([str(python), "-m", "pip", "show", "biopython"],
                           capture_output=True, text=True)
    if shown.returncode != 0:
        subprocess.run([str(python), "-m", "pip", "install", "-q",
                        "--disable-pip-version-check", "biopython"], check=True)
        shown = subprocess.run([str(python), "-m", "pip", "show", "biopython"],
                               capture_output=True, text=True, check=True)
    for line in shown.stdout.splitlines():
        if line.startswith("Version:"):
            return line.split(":", 1)[1].strip()
    return "unknown"


def cached_genomes(count: int) -> list[tuple[str, str]]:
    """Genomes from the fetch cache, so this needs no network.

    Deterministic: the files are sorted and sampled with a fixed seed, so two runs
    compare the same genomes.
    """
    files = sorted(CACHE.glob("*.gz"))
    if not files:
        raise SystemExit(f"no cached records in {CACHE} - run scripts/fetch.py first")
    rng = random.Random(SEED)
    out: list[tuple[str, str]] = []
    for path in rng.sample(files, min(len(files), count * 3)):
        with gzip.open(path, "rb") as fh:
            text = fh.read().decode("utf-8", "replace")
        sequence = "".join(
            c for c in text.split("ORIGIN", 1)[-1].upper() if c in "ACGTRYSWKMBDHVN"
        ) if "ORIGIN" in text else ""
        if len(sequence) > 25_000:
            out.append((path.stem, sequence))
        if len(out) == count:
            break
    if not out:
        raise SystemExit("no cached record held a sequence over 25,000 bases")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genomes", type=int, default=20)
    ap.add_argument("--out", type=Path, default=ROOT / "docs" / "aligner-comparison.json")
    ap.add_argument("--refresh", action="store_true",
                    help="re-run the Biopython side even if its cached result matches "
                         "this run's genomes, seed and biopython version")
    args = ap.parse_args()

    version = ensure_biopython()
    genomes = cached_genomes(args.genomes)
    oligos = [(assay.name, o) for assay in ASSAYS for o in assay.oligos]
    print(f"biopython {version}: {len(oligos)} oligo(s) x {len(genomes)} genome(s) "
          f"= {len(oligos) * len(genomes)} pair(s)\n", flush=True)

    work, mine = [], {}
    start = time.perf_counter()
    for accession, sequence in genomes:
        for assay_name, oligo in oligos:
            hit = find_oligo(sequence, oligo.sequence, oligo.role)
            key = f"{accession}|{oligo.name}"
            searched_for = (
                reverse_complement(oligo.sequence)
                if oligo.role == "reverse"
                else oligo.sequence
            )
            # The same site, counted two ways: with IUPAC codes read as sets (what this
            # tool does) and with plain letter equality (what a substitution matrix of
            # +1/-1 does). One variable between them.
            plain = None
            if hit.found:
                window = sequence[hit.position : hit.position + len(searched_for)]
                if len(window) == len(searched_for):
                    plain = sum(1 for a, b in zip(searched_for, window) if a != b)
            mine[key] = {
                "assay": assay_name,
                "oligo": oligo.name,
                "role": oligo.role,
                "accession": accession,
                "found": hit.found,
                "position": hit.position if hit.found else None,
                "mismatches": hit.mismatches if hit.found else None,
                "ambiguous": hit.ambiguous if hit.found else None,
                "unknown": getattr(hit, "unknown", None),
                "mismatches_plain_identity": plain,
                "degenerate_bases_in_oligo": sorted(set(oligo.sequence) - set("ACGT")),
            }
            # The same string both tools search for. `find_oligo` reverse-complements a
            # reverse primer, because that is the strand a forward-strand genome record
            # carries; handing Biopython the oligo as written asked it a different
            # question and made all five reverse primers look like disagreements.
            searched = (
                reverse_complement(oligo.sequence)
                if oligo.role == "reverse"
                else oligo.sequence
            )
            work.append({"key": key, "oligo": searched, "accession": accession})
    ours = time.perf_counter() - start
    print(f"assay-drift: {len(work)} pair(s) in {ours:.1f}s", flush=True)

    # The expensive half is deterministic given (genomes, seed, biopython version), and
    # it takes minutes. It is cached beside the artefact so that iterating on the REPORT
    # is free - which is the only reason a KeyError in the formatting cost 232 seconds of
    # alignment once.
    cache = args.out.with_suffix(".biopython-cache.json")
    stamp = {"version": version, "seed": SEED,
             "keys": sorted(item["key"] for item in work)}
    rows = None
    if cache.exists() and not args.refresh:
        try:
            held = json.loads(cache.read_text(encoding="utf-8"))
            if held.get("stamp") == stamp:
                rows = held["rows"]
                print(f"biopython:   {len(rows)} pair(s) from {cache.name}", flush=True)
        except (OSError, ValueError, KeyError):
            rows = None

    theirs = 0.0
    start = time.perf_counter()
    with tempfile.TemporaryDirectory() as tmp:
        payload = Path(tmp) / "pairs.json"
        if rows is None:
            payload.write_text(
                json.dumps({"targets": dict(genomes), "pairs": work}), encoding="utf-8"
            )
            proc = subprocess.run(
                [str(python_in(VENV)), "-c", BIOPYTHON_SIDE, str(payload)],
                capture_output=True, text=True,
            )
            theirs = time.perf_counter() - start
            marker = proc.stdout.find("__JSON__")
            if marker < 0:
                raise SystemExit(
                    f"biopython side failed:\n{(proc.stderr or proc.stdout)[-800:]}"
                )
            rows = json.loads(proc.stdout[marker + len("__JSON__"):].splitlines()[0])
            print(f"biopython:   {len(rows)} pair(s) in {theirs:.1f}s", flush=True)
            cache.write_text(json.dumps({"stamp": stamp, "seconds": round(theirs, 1),
                                         "rows": rows}, indent=1), encoding="utf-8")
        else:
            theirs = float(
                json.loads(cache.read_text(encoding="utf-8")).get("seconds") or 0.0
            )

    agree_site, agree_mm, only_mine, only_theirs, differ = [], [], [], [], []
    partial_rows: list[str] = []
    for row in rows:
        key = row["key"]
        ours_row = mine[key]
        if "error" in row:
            only_mine.append({**ours_row, "biopython": row["error"]})
            continue
        if not ours_row["found"]:
            only_theirs.append({**ours_row, "biopython_position": row["position"],
                                "biopython_mismatches": row["mismatches"]})
            continue
        partial = row["unaligned_oligo_bases"] > 0
        same_site = ours_row["position"] == row["position"]
        if same_site:
            agree_site.append(key)
            if ours_row["mismatches"] == row["mismatches"]:
                agree_mm.append(key)
            else:
                differ.append({**ours_row, "biopython_mismatches": row["mismatches"],
                               "biopython_unaligned": row["unaligned_oligo_bases"]})
        else:
            differ.append({**ours_row, "biopython_position": row["position"],
                           "biopython_mismatches": row["mismatches"],
                           "biopython_unaligned": row["unaligned_oligo_bases"]})
        if partial:
            partial_rows.append(row["key"])

    out = {
        "measured": time.strftime("%Y-%m-%d"),
        "biopython_version": version,
        "seed": SEED,
        "genomes": len(genomes),
        "oligos": len(oligos),
        "pairs": len(rows),
        "assay_drift_seconds": round(ours, 1),
        "biopython_seconds": round(theirs, 1),
        "both_found_the_same_site": len(agree_site),
        "and_agreed_on_the_mismatch_count": len(agree_mm),
        "disagreed": len(differ),
        "biopython_left_oligo_bases_unaligned": len(partial_rows),
        "found_only_by_biopython": len(only_theirs),
        "biopython_errors": len(only_mine),
        "iupac_explains": sum(
            1 for r in differ
            if r.get("mismatches_plain_identity") is not None
            and r.get("biopython_mismatches") == r["mismatches_plain_identity"]
            and r["mismatches_plain_identity"] != r["mismatches"]
        ),
        "disagreements": differ[:200],
        "found_only_by_biopython_rows": only_theirs[:200],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1), encoding="utf-8")

    print(f"\n  same site                  {len(agree_site):>6} of {len(rows)}")
    print(f"    and same mismatch count  {len(agree_mm):>6}")
    print(f"  disagreed                  {len(differ):>6}")
    print(f"  found only by biopython    {len(only_theirs):>6}")
    print(f"  biopython aligned only part of the oligo  {len(partial_rows):>6}")
    explained = out["iupac_explains"]
    print(f"  ... of which biopython's count equals plain letter equality at the same "
          f"site: {explained}")
    for row in differ[:8]:
        print(f"    {row['accession'][:12]} {row['oligo']:<8} "
              f"mine {row['mismatches']}mm at {row['position']} | "
              f"biopython {row.get('biopython_mismatches')}mm at "
              f"{row.get('biopython_position', row['position'])}"
              f" ({row.get('biopython_unaligned')} oligo base(s) unaligned)")
    print(f"\nwritten to {args.out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
