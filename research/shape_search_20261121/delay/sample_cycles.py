"""Sample the certified delay sinks in the unscaled original delay coordinates."""
import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from period_windows import parameters
from sink import step
from search import save


def collect(source,out):
    paths=["period-39.json","period-windows.json","high-period-windows.json","long-cycle-blocks.json"]
    mp.mp.dps=90
    out.mkdir(parents=True,exist_ok=True)
    results=[]
    for name in paths:
        data=json.loads((source/name).read_text())
        if isinstance(data,dict):
            data=[data]
        for v in data:
            period=v["period"]
            if "parameters" in v:
                p=list(map(mp.mpf,v["parameters"]))
            else:
                p=parameters(mp.mpf(v["t_numerator"])/v["t_denominator"])
            x=list(map(mp.mpf,v["center"]))
            orbit=[]
            for _ in range(period):
                orbit.append(list(map(float,x)))
                x,_=step(x,p,mp.exp)
            error=float(max(abs(x[i]-mp.mpf(v["center"][i])) for i in range(3)))
            assert error<1e-62,error
            one=np.asarray(orbit,dtype=float)
            # Each sample is the same full periodic geometry, traversed 16 times.
            # This preserves equal phase weights in each temporal half.
            cloud=np.tile(one,(16,1))
            path=out/f"period-{period}.npz"
            np.savez_compressed(path,points=cloud)
            results.append({"period":period,"certificate":name,"traversals":16,
                            "equally_weighted_phases":period,"shape":list(cloud.shape),
                            "closure_at_90_decimal_digits":error,
                            "file_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
                            "cloud_relative":"delay/sink-clouds/"+path.name})
    save(source/"sink-clouds.json",results)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--certificates",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    collect(args.certificates,args.out)
