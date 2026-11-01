"""Archive original calculations and supply relocatable reproduction scripts."""
from pathlib import Path
import shutil,json,platform,importlib.metadata,hashlib

scratch=Path(__file__).parent;repo=Path('D:/Codex/chaos-attractor-lab');out=repo/'results/unexpected-search';code=repo/'research/unexpected_search'
out.mkdir(exist_ok=True);code.mkdir(exist_ok=True);original=out/'source-original';original.mkdir(exist_ok=True)
scripts=['anderson_pilot','anderson_extended','anderson_gallery','anderson_long','anderson_spectra','anderson_precision','anderson_basins','tetra_vortex','tetra_gallery','tetra_long','long_galleries','tetra_spectra','tetra_dop853','tetra_basins','tetra_render','final_visuals','symmetry_audit','compare_unexpected','summarize_evidence']
for name in scripts:
    src=scratch/f'{name}.py';shutil.copy2(src,original/src.name);s=src.read_text()
    s=s.replace('Path(__file__).parent','OUTPUT')
    s=s.replace("str(Path('D:/Codex/chaos-attractor-lab/research'))",'str(RESEARCH)')
    s=s.replace("'D:/Codex/chaos-attractor-lab/research'",'str(RESEARCH)')
    s=s.replace("Path('D:/Codex/chaos-attractor-lab')",'ROOT')
    (code/src.name).write_text('from paths import OUTPUT, ROOT, RESEARCH\n'+s,encoding='utf-8')
for folder in ['anderson-v1','anderson-v2','anderson-long','anderson-spectra','anderson-precision','anderson-basins','tetra-v1','tetra-long','tetra-spectra','tetra-dop853','tetra-basins','tetra-render','comparison','visuals']:
    shutil.copytree(scratch/folder,out/folder,dirs_exist_ok=True,ignore=shutil.ignore_patterns('dense-*.npy'))
shutil.copy2(scratch/'facts.json',out/'facts.json')
shutil.copy2(scratch/'archive_packet.py',original/'archive_packet.py')
for name in ['twist_map.py','render_twist.py','compare_attractor_shapes.py','compare_affine_shapes.py']:
    shutil.copy2(repo/'research'/name,original/name)
for name in ['tetrahedral-flow.png','tetrahedral-flow-oblique.png','tetrahedral-vortex-565.png','tetrahedral-flow.gif','anderson-memory.png']:
    shutil.copy2(scratch/'visuals'/name,repo/'docs/assets'/name)
environment={'python':platform.python_version(),'platform':platform.platform(),'packages':{n:importlib.metadata.version(n) for n in ['numpy','scipy','numba','matplotlib','pillow','mpmath']},'original_scratch':str(scratch),'repository_baseline':'6bbaded','notes':['All new calculation and temporary files were on drive D.','Portable scripts change only path resolution; exact original source snapshots are included.','The first dense rendering attempt failed after numerical integration due to unpacking the colorize return tuple incorrectly. The corrected script reran the integration, saved dynamics before rendering, and produced the archived images.','The first interval-root attempt rejected NumPy scalar repr strings before any interval claim or basin result was produced. Scalar conversion was corrected and the complete run succeeded.','The dense arrays omitted from Git are reproducible; rendering histograms, hashes, sample clouds and all dynamics summaries are retained.']}
(out/'environment.json').write_text(json.dumps(environment,indent=2))
reference_files=sorted((repo/'results/attractor-novelty/reference/trajectories').glob('*.npz'))+[repo/'results/twist-geometry/long-runs.npz',repo/'results/gaussian-vortex/long-runs.npz']
(out/'reference-inputs.json').write_text(json.dumps([{'path':p.relative_to(repo).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in reference_files],indent=2))
print(json.dumps({'portable_scripts':len(scripts),'archived_bytes':sum(p.stat().st_size for p in out.rglob('*') if p.is_file())}),flush=True)
