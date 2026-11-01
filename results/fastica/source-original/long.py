from pathlib import Path
import json,sys,time,numpy as np
root=Path('D:/Codex/unexpected-attractors-20261101/fastica');sys.path.insert(0,str(root));from follow import integrate
out=root/'long-01';out.mkdir(exist_ok=False);rows=[]
for case in [15,40,204,246,351]:
 data=json.loads((root/f'follow-01/case-{case}-input.json').read_text());X=np.array(data['X']);w=np.array(data['start']);start=time.time()
 rates,cloud,mn,mg,last=integrate(X,w,data['kind'],data['alpha'],10000,5000000,500000,50)
 np.savez_compressed(out/f'case-{case}.npz',rates=rates,cloud=cloud,endpoint=last)
 row={'case':case,'warmup':10000,'observed':5000000,'window':500000,'stride':50,'rates':rates.tolist(),'mean':rates.mean(axis=0).tolist(),'minimum_update_norm':mn,'minimum_qr_diagonal':mg,'sphere_error':float(abs(np.sum(cloud*cloud,axis=1)-1).max()),'wall_seconds':time.time()-start};rows.append(row);print(json.dumps(row),flush=True)
 (out/'report.json').write_text(json.dumps(rows,indent=2)+'\n')
