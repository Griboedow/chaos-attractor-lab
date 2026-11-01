"""Declared time-step refinement, with a new IVP startup before tangent burn."""
from pathlib import Path
import subprocess
import sys

root=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(root/'restart_accelerated.py'),
    '--seed',str(root/'spectrum124-n64'),'--output',str(root/'field124-dt003-relax'),
    '--segments','64','--dt','0.003','--duration','24'],check=True,cwd=root)
subprocess.run([sys.executable,str(root/'tangent.py'),
    '--seed',str(root/'field124-dt003-relax'),'--output',str(root/'spectrum124-dt003'),
    '--duration','120','--burn','24'],check=True,cwd=root)
