#!/usr/bin/env python3
"""Compare all real-data tables with the audited pre-refactor baseline."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
expected = json.loads((root / "tests/expected_tables.json").read_text())
actual = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
          for p in (root / "outputs/tables").glob("*.csv")}
if actual != expected:
    changed = sorted(k for k in set(actual) | set(expected) if actual.get(k) != expected.get(k))
    raise SystemExit("Tables differ from audited baseline: " + ", ".join(changed))
print(f"All {len(expected)} tables reproduce the audited baseline exactly.")
