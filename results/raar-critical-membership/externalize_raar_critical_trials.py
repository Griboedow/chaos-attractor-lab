"""Preserve large proposal-only arrays in the external cache before publication."""
import json
import shutil
from pathlib import Path
import numpy as np
from audit_raar_coverings import sha


def main():
    repo=Path(__file__).resolve().parents[1];packet=(repo/'results/raar-critical-membership').resolve();cache=Path('D:/Codex/chaos-attractor-cache')
    retained='runs/raar-critical-membership-render02/figure-data.npz';rows=[]
    for source in sorted(packet.rglob('*.npz')):
        relative=source.relative_to(packet).as_posix()
        if source.stat().st_size<262144 or relative==retained:continue
        assert source.resolve().is_relative_to(packet)
        digest=sha(source);destination=cache/digest[:2]/(digest+'.npz');destination.parent.mkdir(parents=True,exist_ok=True)
        if not destination.exists():shutil.copyfile(source,destination)
        assert sha(destination)==digest and destination.stat().st_size==source.stat().st_size
        with np.load(source,allow_pickle=False) as data:
            arrays=[{'name':name,'dtype':str(data[name].dtype),'shape':list(data[name].shape),'bytes':data[name].nbytes,'retained':False} for name in data.files]
        command=('python research/propose_raar_thin_window.py --out <new-directory>' if 'thin-window' in relative else 'python research/propose_raar_critical_crossing.py --out <new-directory>')
        row={'path':source.relative_to(repo).as_posix(),'bytes':source.stat().st_size,'sha256':digest,'git_blob':None,'group':'raar-critical-membership','compact':None,'arrays':arrays,'local_cache':str(destination),'recovery':'recompute','recompute_command':command,'scope':'Optional numerical trial array; not required by the fold-membership proof.'}
        if 'render01' in relative:
            row.update(recovery='retained_duplicate',retained_source='results/raar-critical-membership/'+retained,recompute_command='python research/render_raar_critical_membership.py --proof results/raar-critical-membership/runs/raar-critical-thin-fold-region01 --out <new-directory>')
        rows.append(row)
        source.unlink()
    (packet/'omissions.json').write_text(json.dumps({'threshold_bytes':262144,'files':rows,'retained_large_input':{'path':'results/raar-critical-membership/'+retained,'bytes':(packet/retained).stat().st_size,'reason':'Exact complete input to the published fold-membership diagnostic and its bitwise figure audit; 8,193 guide states, determinants, camera and certified location hull.'}},indent=2)+'\n',encoding='utf-8')
    (packet/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print(json.dumps({'omitted_files':len(rows),'external_bytes':sum(r['bytes'] for r in rows),'retained_figure_bytes':(packet/retained).stat().st_size}),flush=True)


if __name__=='__main__':main()
