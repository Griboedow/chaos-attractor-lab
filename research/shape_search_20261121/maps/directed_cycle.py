"""Prove a strict cyclic return using separate directed-interval partitions."""
import argparse
import hashlib
import json
from pathlib import Path
from directed_cover import region


def verify(packet,path,budget,parameter_bits=None):
    period=packet["case_order"]
    checks=[(step,0,"avoid") for step in range(1,period)]
    checks.extend((period,j,"coordinate") for j in range(3))
    digest=hashlib.sha256(json.dumps(packet,sort_keys=True).encode()).hexdigest()
    results=[]
    for i,(step,axis,kind) in enumerate(checks):
        partial=path.with_name(path.stem+f"-part-{i}.json")
        if partial.exists():
            saved=json.loads(partial.read_text())
            if saved["box_digest"]!=digest or saved.get("parameter_radius_bits")!=parameter_bits:
                raise RuntimeError("Incompatible partial result")
            info=saved["check"]
        else:
            cover=dict(axis=axis,faces=[-1.,1.],map_steps=step)
            info=region(packet,cover,[-1.]*3,[1.]*3,kind,None,None,budget,parameter_bits)
            info.update(map_steps=step,axis=axis)
            partial.write_text(json.dumps(dict(box_digest=digest,parameter_radius_bits=parameter_bits,check=info))+"\n")
        results.append(info)
    proof=dict(case=packet["case"],phase_regions=period,box_digest=digest,parameter_radius_bits=parameter_bits,
               method="separate coordinate and phase-exclusion checks with directed binary64 intervals",checks=results)
    path.write_text(json.dumps(proof,indent=2)+"\n")
    print(packet["case"],period,[r["nodes"] for r in results],flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("box",type=Path)
    ap.add_argument("--output",type=Path,required=True);ap.add_argument("--budget",type=int,default=400_000)
    ap.add_argument("--parameter-radius-bits",type=int)
    a=ap.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
    verify(json.loads(a.box.read_text()),a.output,a.budget,a.parameter_radius_bits)
