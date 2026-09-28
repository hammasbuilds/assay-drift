"""The 2019 reference genome, vendored so nothing here needs the network.

NC_045512.2 (Wuhan-Hu-1) is the genome every assay in `primers.py` was designed
against, so it is what "drift" is measured *from*: an oligo's position in this
sequence is the coordinate a mutation gets named by, and an oligo that does not
match it is a typo rather than a discovery.

It ships inside the package (9 KB gzipped) rather than beside the tests,
because three different callers need it - the test suite, `demo.py`, and the
`check` / `mutations` commands - and only one of those still exists after
`pip install`. Its sha256 is pinned below: a corrupted or swapped fixture would
silently change what "matches the reference" means, which is the one failure
this project could not detect from its own output.

The reference has not changed since 2020, so vendoring it costs nothing in
freshness and buys a suite that passes on a plane, behind a proxy, and on the
day NCBI has an outage.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from functools import lru_cache
from pathlib import Path

ACCESSION = "NC_045512.2"
PATH = Path(__file__).resolve().parent / "data" / "NC_045512.2.json.gz"

# Pinned so a corrupted or swapped fixture fails loudly.
SHA256 = "7d5621cd3b3e498d0c27fcca9d3d3c5168c7f3d3f9776f3005c7011bd90068ca"
LENGTH = 29903


class ReferenceError(RuntimeError):
    """The vendored reference genome is missing or is not the expected one."""


@lru_cache(maxsize=1)
def genome() -> str:
    """The reference sequence, verified against its pinned digest."""
    if not PATH.exists():
        raise ReferenceError(f"vendored reference genome missing at {PATH}")
    with gzip.open(PATH, "rt", encoding="utf-8") as handle:
        sequence = json.load(handle)["sequence"]
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    if digest != SHA256:
        raise ReferenceError(
            f"{PATH} is not {ACCESSION}: sha256 {digest}, expected {SHA256}. "
            f"Re-fetch it rather than updating the digest - the reference has "
            f"not changed since 2020."
        )
    return sequence
