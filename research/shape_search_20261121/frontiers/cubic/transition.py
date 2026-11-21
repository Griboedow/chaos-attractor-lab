"""Survey basin fates through the knotting transition of a primitive-period cubic drive."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist
from gradient import independent,integrate
from model import forcing,geometric_word,invariants,xyz


def distinct_sets(base,ends):
    records=[]
    for i in range(ends.shape[1]):
        p=next((k for k in range(1,17) if np.max(abs(ends[-5:,i]-ends[-5-k:-k,i]))<2e-7),None)
        if p is None:continue
        phase=ends[-p:,i]
        for rec in records:
            if rec["period"]!=p:continue
            other=ends[-p:,rec["representative"]]
            nearest=abs(phase[:,None]-other[None,:])
            if max(nearest.min(0).max(),nearest.min(1).max())<5e-6:
                rec["starts"].append(i);break
        else: records.append(dict(period=p,representative=i,starts=[i]))
    return records


def survey(row,out,packet,steps,turns):
    rng=np.random.default_rng(4521)
    starts=np.r_[[-np.sqrt(3),0.,np.sqrt(3)],rng.uniform(-2,2,12)+1j*rng.uniform(-1.6,1.6,12)]
    levels=sorted(set([.03,.15,.23,.34]+[round(x,3) for x in np.arange(.075,.351,.005)]))
    b=forcing(np.linspace(0,2*np.pi,2*steps+1),row)
    result=[];receipts=[]
    for gain in levels:
        paths,ends=integrate(starts,b,gain,turns)
        distinct=distinct_sets(starts,ends)
        entry=dict(gain=gain,components=distinct,
                   unclassified=len(starts)-sum(len(q["starts"]) for q in distinct))
        reps=[q["representative"] for q in distinct if q["period"]==1]
        if len(reps)==3:
            sample=paths[:,reps]
            values=[invariants(geometric_word(sample,x)) for x in (.181,.576)]
            entry["geometric_views"]=values
        if gain in (.03,.15,.23,.34):
            controls=[]
            for rec in distinct:
                k=rec["period"];j=rec["representative"]
                phase=ends[-1,j];segments=[];monodromy=np.eye(2)
                start=phase
                for _ in range(k):
                    curve,m=independent(row,gain,[phase.real,phase.imag],full_matrix=True)
                    segments.append(curve);phase=curve[-1];monodromy=m@monodromy
                controls.append(dict(period=k,closure=float(abs(phase-start)),
                                     least_phase_distance=min([float(abs(z[0]-start)) for z in segments[1:]] or [0.]),
                                     dominant_multiplier_modulus=float(max(abs(np.linalg.eigvals(monodromy))))))
                offsets=np.linspace(0,2*np.pi,4097)
                for piece in segments:receipts.append((gain,j,xyz(offsets,piece,row)))
            entry["independent_DOP853"]=controls
        result.append(entry)
    out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out/"transition-curves.npz",points=np.vstack([v for _,_,v in receipts]),
                        offsets=np.cumsum([0]+[len(v) for _,_,v in receipts]),
                        gain=[g for g,_,_ in receipts],orbit=[j for _,j,_ in receipts])
    from hashlib import sha256
    data=dict(case=row["case"],turns=turns,steps_per_turn=steps,initials=np.column_stack((starts.real,starts.imag)).tolist(),
              control_cloud="frontiers/cubic/transition-curves.npz",
              sha256=sha256((out/"transition-curves.npz").read_bytes()).hexdigest(),sweep=result)
    packet.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print([(m["gain"],[(r["period"],len(r["starts"])) for r in m["components"]]) for m in result if "independent_DOP853" in m],flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("certificates",type=Path)
    ap.add_argument("--out",type=Path,required=True);ap.add_argument("--manifest",type=Path,required=True)
    ap.add_argument("--steps",type=int,default=2048);ap.add_argument("--turns",type=int,default=64)
    args=ap.parse_args();row=next(r for r in json.loads(args.certificates.read_text())["rows"] if r["case"]==48000001)
    survey(row,args.out,args.manifest,args.steps,args.turns)
