"""Frozen inputs and explicitly requested output for the return-map follow-up."""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parents[2]
INPUT=Path(os.environ.get('CHAOS_ICOSA_INPUT',ROOT/'results/icosahedral-flow'))
STRUCTURE=Path(os.environ.get('CHAOS_ICOSA_STRUCTURE_INPUT',ROOT/'results/icosahedral-structure'))
CYCLES=Path(os.environ.get('CHAOS_ICOSA_CYCLES_INPUT',ROOT/'results/icosahedral-cycles'))
if 'CHAOS_ICOSA_RETURNS_OUTPUT' not in os.environ:
    raise RuntimeError('Set CHAOS_ICOSA_RETURNS_OUTPUT to a scratch directory; stages refuse existing run folders.')
OUTPUT=Path(os.environ['CHAOS_ICOSA_RETURNS_OUTPUT'])
OUTPUT.mkdir(parents=True,exist_ok=True)
