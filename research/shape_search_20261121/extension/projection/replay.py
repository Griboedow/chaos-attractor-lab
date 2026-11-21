"""Replay a retained configuration or extend its QR windows."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from confirm import from_record
from affine_model import orbit

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--record',type=Path,required=True);p.add_argument('--cloud',type=Path);p.add_argument('--extended',type=Path);p.add_argument('--reference',type=Path)
    a=p.parse_args();rec=json.loads(a.record.read_text());raw=from_record(rec);initial=np.array(rec['initial'])
    if a.extended:
        tail,qr=orbit(raw,initial,100000,5000000,20,4096)
        a.extended.parent.mkdir(parents=True,exist_ok=True)
        a.extended.write_text(json.dumps(dict(case=rec['case'],burn=100000,window=5000000,windows=qr.tolist(),rates=qr.mean(0).tolist(),minimum=qr.min(0).tolist()),indent=2)+'\n',encoding='utf8')
    else:
        if a.cloud is None:p.error('--cloud or --extended is required')
        z,qr=orbit(raw,initial,rec['burn'],rec['window_length'],rec['windows'],rec['points_saved'])
        a.cloud.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(a.cloud,points=z.astype('f4'),case=rec['case'],initial=initial)
        check=dict(case=rec['case'],max_qr_difference=float(np.max(abs(qr-np.asarray(rec['qr_windows'])))),sha256=hashlib.sha256(a.cloud.read_bytes()).hexdigest())
        if a.reference:check['identical_cloud']=np.array_equal(z.astype('f4'),np.load(a.reference)['points'])
        print(check)
