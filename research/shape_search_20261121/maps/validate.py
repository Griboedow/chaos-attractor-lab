"""Long trajectories, basin samples, parameter controls and independent spectra."""
import argparse
import hashlib
import json
import sys
from pathlib import Path
from search import CACHE, RESULTS, REPO, FAMILY, START, sample, advance
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0,str(REPO/"research"))
from twist_map import spectra, step, step_tangent


def write_json(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    text=json.dumps(value,indent=2,allow_nan=False)+"\n"
    temp=path.with_name(path.name+".tmp")
    temp.write_text(text,encoding="utf-8")
    temp.replace(path)


def outcome(tree, diagonal, rates, z):
    distances=tree.query(z,workers=1)[0]/diagonal
    return dict(leading=float(rates.mean()),window_range=[float(rates.min()),float(rates.max())],
                nearest_reference_diagonal_quantiles=np.quantile(distances,[.5,.95,.99]).tolist())


def validate(case,label,p,cache,results):
    rng=np.random.default_rng(np.random.SeedSequence([20261121,case,61]))
    r,tail=sample(p[None,:],START[None,:],100_000,2_000_000,40,240_000)
    z=tail[0]
    tree=cKDTree(z)
    diag=np.linalg.norm(np.ptp(z,axis=0))
    if diag<1e-11:
        row=dict(case=case,label=label,family=FAMILY,params=p.tolist(),primary_start=START.tolist(),
                 excluded="Reference orbit collapsed after long burn-in.",reference_rates=r[0].tolist())
        write_json(results/"checks"/f"{case}.json",row)
        return
    starts=np.vstack((START,rng.uniform(-2,2,(15,3))))
    para=np.tile(p,(len(starts),1))
    w,b=sample(para,starts,60_000,250_000,10,16_384)
    basin=[dict(start=s.tolist(),**outcome(tree,diag,ww,zz)) for s,ww,zz in zip(starts,w,b)]
    offsets=rng.uniform(-1e-3,1e-3,(8,10))
    offsets[:,7:10][:,p[7:10]==0]=0
    near=np.vstack((p+offsets,np.round(p,3)))
    wc,bc=sample(near,np.tile(START,(len(near),1)),50_000,200_000,8,8192)
    cohort=[dict(params=q.tolist(),**outcome(tree,diag,ww,zz)) for q,ww,zz in zip(near,wc,bc)]
    # Reuse the original model's independently authored stage implementation,
    # full tangent, and twice-reorthogonalized QR. The scout uses its own Jv.
    windows,ind,_,_,radii=spectra(p[None,:],START[None,:],80_000,600_000,30_000,32_768)
    trace=3*np.log(p[0])
    algebra=[]
    for trial in rng.uniform(-1,1,(12,3)):
        native, jac=step_tangent(trial.copy(),np.eye(3),p)
        v=rng.normal(size=3)
        v/=np.linalg.norm(v)
        ztest,vtest=trial.copy(),v.copy()
        growth=advance(p,ztest,vtest)
        algebra.append([float(np.max(abs(ztest-native))),float(np.max(abs(vtest*np.exp(growth)-jac@v))),
                        float(abs(np.linalg.det(jac)-p[0]**3))])
    pts=z.astype(np.float32)
    cloud=cache/"clouds"/f"{case}.npz"
    cloud.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(cloud,points=pts)
    row=dict(case=case,label=label,family=FAMILY,parameter_order=["r","a_xy","a_yz","a_zx","k_xy","k_yz","k_zx","h_xy","h_yz","h_zx"],
        params=p.tolist(),primary_start=START.tolist(),reference_rates=r[0].tolist(),
        reference_mean=float(r.mean()),independent_spectrum=windows[0].mean(axis=0).tolist(),
        independent_spectrum_windows=windows[0].tolist(),maximum_trace_error=float(np.max(abs(windows[0].sum(axis=1)-trace))),
        max_state_tangent_determinant_errors=np.max(algebra,axis=0).tolist(),
        independent_tail=outcome(tree,diag,windows[0,:,0],ind[0]),native_bounding_box_diagonal=float(diag),
        maximum_reference_radius=float(np.linalg.norm(z,axis=1).max()),maximum_independent_radius=float(radii[0]),
        basin=basin,local_parameters=cohort,cloud=str(cloud.resolve()),key="points",dtype="float32",points=len(z),
        sha256=hashlib.sha256(cloud.read_bytes()).hexdigest())
    write_json(results/"checks"/f"{case}.json",row)
    print(case,label,row["independent_spectrum"],sum(x["window_range"][0]>0 for x in basin),
          sum(x["window_range"][0]>0 for x in cohort),flush=True)


if __name__=="__main__":
    args=argparse.ArgumentParser(description=__doc__)
    args.add_argument("--cache",type=Path,default=CACHE)
    args.add_argument("--results",type=Path,default=RESULTS)
    args.add_argument("--preview",type=Path,required=True)
    args.add_argument("--proposals",type=Path,required=True)
    args.add_argument("--cases",type=int,nargs="*")
    a=args.parse_args()
    lookup={x["case"]:np.array(x["params"]) for x in json.loads(a.preview.read_text())}
    proposed=json.loads(a.proposals.read_text())
    for case,label in proposed.items():
        case=int(case)
        if a.cases and case not in a.cases:
            continue
        validate(case,label,lookup[case],a.cache,a.results)
