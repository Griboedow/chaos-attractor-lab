"""Resolve the repository and require a separate output location for reproduction."""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parents[2]
RESEARCH=ROOT/'research'
OUTPUT=Path(os.environ.get('CHAOS_ICOSA_OUTPUT', ROOT/'results/icosahedral-flow'))
OUTPUT.mkdir(parents=True,exist_ok=True)
