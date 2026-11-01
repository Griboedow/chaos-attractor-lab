"""Make full final-window portrait inputs by translation and proper rotation only."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);source=args.input/'3x2-final-250000-states.npz';a=np.load(source);records=[]
    e,az=np.deg2rad([27,-45]);projection=np.array([[-np.sin(az),np.cos(az),0],[-np.cos(az)*np.sin(e),-np.sin(az)*np.sin(e),np.cos(e)]])
    for case,z in zip(a['case'],a['states']):
        assert z.shape==(250000,3) and np.isfinite(z).all();center=z.mean(0);u,s,vh=np.linalg.svd(z-center,full_matrices=False);basis=vh.T
        if np.linalg.det(basis)<0:basis[:,-1]*=-1
        points=(z-center)@basis;path=args.output/f'case-{case}-points.npz';np.savez_compressed(path,points=points)
        record=dict(case=int(case),samples=len(z),start_step=105790001,stop_step=106040000,native_source=str(source),native_source_sha256=digest(source),native_center=center.tolist(),native_basis=basis.tolist(),pca_variance=(s*s/np.sum(s*s)).tolist(),basis=np.eye(3).tolist(),center=((points.min(0)+points.max(0))/2).tolist(),projection=projection.tolist(),point_file_sha256=digest(path),geometry='All 250000 consecutive states; translation and proper orthogonal PCA rotation, without rescaling, smoothing, replication, or connecting discrete iterates.')
        assert np.max(abs(basis.T@basis-np.eye(3)))<1e-12
        (args.output/f'case-{case}-view.json').write_text(json.dumps(record,indent=2)+'\n');records.append(record)
    print(json.dumps([dict(case=r['case'],variance=r['pca_variance']) for r in records]))
