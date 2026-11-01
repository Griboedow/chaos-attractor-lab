"""Explicit D-drive-compatible output, with separately selectable old inputs."""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parents[2]
RESEARCH=ROOT/'research'
INPUT=Path(os.environ.get('CHAOS_ICOSA_INPUT', ROOT/'results/icosahedral-flow'))
if 'CHAOS_ICOSA_STRUCTURE_OUTPUT' not in os.environ:
    raise RuntimeError('Set CHAOS_ICOSA_STRUCTURE_OUTPUT to a new scratch directory')
OUTPUT=Path(os.environ['CHAOS_ICOSA_STRUCTURE_OUTPUT'])
OUTPUT.mkdir(parents=True,exist_ok=True)
