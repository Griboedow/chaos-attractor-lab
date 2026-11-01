"""Compare every completed common chunk before and after early subdivision."""
import argparse
import json
from pathlib import Path
import numpy as np
from audit_raar_trapping_batch import file_hash as sha


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    parts=[]
    for i in range(4):
        old=args.base/f'raar-trapping-audit-round1-part{i}'
        new=args.base/f'raar-trapping-audit-round2-part{i}'
        assert (old/'report.json').exists() or (old/'interrupted.json').exists()
        checkpoint=old/('completed-checkpoint.json' if (old/'interrupted.json').exists() else 'progress.json')
        a=json.loads(checkpoint.read_text());b=json.loads((new/'report.json').read_text())
        declaration=json.loads((old/'declaration.json').read_text())
        for name,digest in declaration['source_sha256'].items():
            assert sha(old/name)==digest and sha(Path(__file__).with_name(name))==digest
        assert b['completed']
        assert len(a)<=len(b['history'])
        for x,y in zip(a,b['history']):
            assert {k:v for k,v in x.items() if k!='seconds'}=={k:v for k,v in y.items() if k!='seconds'}
            assert sha(old/x['new_file'])==sha(new/y['new_file'])
        parts.append({'part':i,'compared_sources':sum(v['processed'] for v in a),'compared_chunks':len(a),
                      'old_completed':(old/'report.json').exists(),'old_checkpoint_file':checkpoint.name,'old_progress_sha256':sha(checkpoint),
                      'new_report_sha256':sha(new/'report.json')})
    a=np.load(args.base/'raar-interval-grid-validation05/sample.npz');b=np.load(args.base/'raar-interval-grid-validation06/sample.npz')
    assert a.files==b.files and all(np.array_equal(a[k],b[k]) for k in a.files)
    result={'status':'identical_common_results','parts':parts,'total_compared_sources':sum(v['compared_sources'] for v in parts),
            'validation_arrays_identical':a.files,'source_sha256':sha(Path(__file__)),
            'scope':'Two raster control flows with the same outward arithmetic and enclosures. This is not an independent arithmetic proof stack. Only completed common chunks are compared; incomplete work is excluded.'}
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
