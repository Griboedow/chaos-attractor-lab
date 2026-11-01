"""Separate frozen inputs from explicit scratch output."""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parents[2]
RESEARCH=ROOT/'research'
INPUT=Path(os.environ.get('CHAOS_ICOSA_INPUT', ROOT/'results/icosahedral-flow'))
STRUCTURE=Path(os.environ.get('CHAOS_ICOSA_STRUCTURE_INPUT', ROOT/'results/icosahedral-structure'))
if 'CHAOS_ICOSA_CYCLES_OUTPUT' not in os.environ:
    raise RuntimeError('Set CHAOS_ICOSA_CYCLES_OUTPUT to a new scratch directory')
OUTPUT=Path(os.environ['CHAOS_ICOSA_CYCLES_OUTPUT'])
OUTPUT.mkdir(parents=True,exist_ok=True)
