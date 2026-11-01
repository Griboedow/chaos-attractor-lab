"""Reconstruct the declared design and audit the higher-order controls."""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from scipy.stats import qmc
from cohort_product_alphabet import cohort


def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main(a):
    a.out.mkdir(parents=True,exist_ok=False);set_num_threads(4);begin=time.time()
    for name in [Path(__file__).name,'cohort_product_alphabet.py','search_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']:
        (a.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(a.scout/'screen.npz') as f:
        params=f['params'];levels=f['levels'];orders=f['orders'];flags=f['flags']
    protocol=read(a.scout/'protocol.json');generated=[];controls=[]
    for old in [988,489]:
        p=a.controls/f'case-{old}-states.npz'
        with np.load(p) as f:
            q=int(f['order']);generated.append((q,f['params'],f['levels'][:,:q]))
        controls.append(dict(earlier_case=old,source_sha256=sha(p)))
    for q,beta,half,ratio in itertools.product([5,7,9],[.90,.95,.99],[.45,.75,1.05],[.03,.15,.75]):
        generated.append((q,np.array([beta,ratio*(2*half/(q-1))**2]),np.tile(half*np.linspace(-1,1,q),(3,1))))
    for u in qmc.Sobol(8,scramble=True,seed=2026110152).random_base2(8):
        beta=.85+.149*u[0];ratio=10**(-2+2.6*u[1]);half=.15+1.1*u[2:5];center=1.6*u[5:8]-.8
        for q in [5,7,9]:
            generated.append((q,np.array([beta,ratio*(2*np.exp(np.log(half).mean())/(q-1))**2]),center[:,None]+half[:,None]*np.linspace(-1,1,q)))
    assert len(generated)==len(params)==851
    for i,(q,p,g) in enumerate(generated):
        assert q==orders[i] and np.array_equal(p,params[i]) and np.array_equal(g,levels[i,:,:q])
        row=protocol['cases'][i]
        assert row['id']==i and row['order']==q and [row['beta'],row['tau']]==p.tolist() and row['levels']==g.tolist()
    assert int(flags.sum())==46
    cohort_report=read(a.cohort/'report.json');cohort_rows=[]
    starts=np.vstack([[.317,-.211,.419],8*qmc.Sobol(3,scramble=True,seed=2026110146).random_base2(5)-4])
    assert np.array_equal(starts,cohort_report['starts'])
    for row in cohort_report['cases']:
        i=row['id'];values=cohort(*params[i],levels[i],orders[i],starts,100000,500000)
        with np.load(a.cohort/f'case-{i}.npz') as f:
            assert np.array_equal(starts,f['starts'])
            for key,v in zip(['rates','tails','status','minima'],values):assert np.array_equal(v,f[key],equal_nan=True),(i,key)
        positive=periodic=fixed=0
        for j,record in enumerate(row['starts']):
            z=values[1][j];scale=max(1.,np.linalg.norm(np.ptp(z,axis=0)))
            periods=[k for k in range(1,513) if np.max(np.linalg.norm(z[k:]-z[:-k],axis=1))<1e-8*scale]
            period=min(periods) if periods else None
            assert period==record['period'] and record['spectrum']==values[0][j].mean(0).tolist()
            assert record['minimum_windows']==values[0][j].min(0).tolist()
            positive+=int(period is None and values[0][j,:,0].min()>0)
            periodic+=int(period is not None);fixed+=int(period==1)
        cohort_rows.append(dict(id=i,positive_nonperiodic=positive,periodic=periodic,fixed=fixed,other=33-positive-periodic,minimum_leading_window=float(values[0][:,:,0].min()),replay='All 33 starts: states, rates, status and radius bitwise exact.'))
        print(json.dumps(cohort_rows[-1]),flush=True)
    precision_rows=[]
    for folder in a.precision:
        r=read(folder/'report.json');assert r['source_sha256']==sha(a.scout/'screen.npz')
        for row in r['cases']:
            i=row['id']
            with np.load(folder/f'case-{i}.npz') as f:
                rates=f['rates'];z=f['states'];assert len(z)==10000
                assert np.array_equal(rates,[float(x) for x in row['windows']])
                assert abs(rates.mean()-row['mean'])<1e-14 and rates.min()==row['minimum_window']
                scale=max(1.,np.linalg.norm(np.ptp(z,axis=0)))
                periods=[k for k in range(1,513) if np.max(np.linalg.norm(z[k:]-z[:-k],axis=1))<1e-8*scale]
                precision_rows.append(dict(id=i,digits=r['digits'],mean=row['mean'],minimum_window=row['minimum_window'],detected_period=min(periods) if periods else None,variance_fractions=row['variance_fractions']))
    catalog=read(a.catalog/'report.json');assert len(catalog['reference_files'])==108 and catalog['exact_affine_control']<1e-12
    for ident,representations in catalog['cases'].items():
        for name,row in representations.items():
            assert row['source_sha256']==sha(a.long/f'case-{ident}.npz')
            assert len(row['all_rankings'])==108-len(row['excluded']) and len(row['nearest'])==3
            assert all(np.isfinite(v['best_local_icp_rms']) for v in row['nearest'])
    report=dict(status='passed',design_reconstructed=851,control_sources=controls,source_protocol_sha256=sha(a.scout/'protocol.json'),cohort_replays=cohort_rows,precision_tail_checks=precision_rows,catalog_cases=len(catalog['cases']),catalog_reference_files=108,seconds=time.time()-begin,scope='Exact same-environment replay and summary accounting. Precision controls are checked for finite returns on their stored binary64 visualization tails, not rerun here. Library alignment is not a global equivalence or priority test.')
    (a.out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['scout','controls','cohort','catalog','long','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--precision',type=Path,nargs='+',required=True);main(p.parse_args())
