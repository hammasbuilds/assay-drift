"""Make `pytest` work from a clone with nothing installed.

The package lives in src/, so `import assaydrift` only resolves after an
install. A project that fetches its own data over the network should not also
demand an install step before its tests run.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
