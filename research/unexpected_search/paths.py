"""Keep reproduction output separate from the frozen evidence packet."""
from pathlib import Path
import os

ROOT = Path(__file__).resolve().parents[2]
RESEARCH = ROOT / 'research'
OUTPUT = Path(os.environ.get('CHAOS_SEARCH_OUTPUT', ROOT / 'results/unexpected-search'))
OUTPUT.mkdir(parents=True, exist_ok=True)
