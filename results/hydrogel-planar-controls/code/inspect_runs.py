"""Sampled proper centerline crossings; not a finite-thickness contact model."""
import json
from pathlib import Path
import numpy as np
from numba import njit

@njit(cache=True)
def crossing_counts(theta):
    n, nt = theta.shape
    counts = np.zeros(nt,dtype=np.int64)
    for k in range(nt):
        x,z = np.zeros(n+1),np.zeros(n+1)
        for i in range(n):
            x[i+1] = x[i]+np.sin(theta[i,k])/n
            z[i+1] = z[i]+np.cos(theta[i,k])/n
        for i in range(n):
            ax,az = x[i+1]-x[i],z[i+1]-z[i]
            for j in range(i+2,n):
                bx,bz = x[j+1]-x[j],z[j+1]-z[j]
                cx,cz = x[j]-x[i],z[j]-z[i]
                det = ax*bz-az*bx
                if abs(det) < 1e-15:
                    continue
                u,v = (cx*bz-cz*bx)/det,(cx*az-cz*ax)/det
                if 1e-10 < u < 1-1e-10 and 1e-10 < v < 1-1e-10:
                    counts[k] += 1
    return counts

def inspect(folder):
    r = json.loads((folder/"report.json").read_text(encoding="utf-8"))
    data = np.load(folder/"trajectory.npz")
    late = data["t"] >= r["duration"]/2
    counts = crossing_counts(np.ascontiguousarray(data["y"][:r["segments"],late]))
    times = data["t"][late]
    return dict(folder=str(folder),parameters=r["parameters"],segments=r["segments"],
                duration=r["duration"],method=r["method"],rtol=r["rtol"],success=r["success"],
                late_samples=len(counts),samples_with_proper_crossing=int(np.sum(counts>0)),
                crossing_fraction=float(np.mean(counts>0)),max_crossings=int(counts.max()),
                first_sampled_late_crossing=float(times[np.flatnonzero(counts>0)[0]]) if np.any(counts) else None,
                period_mean=r["period_mean"],period_std=r["period_std"],
                full_state_return_max_last10=r.get("full_state_return_max_last10"))

if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    paths = (sorted((root/"scout-01").glob("case-*"))+sorted(root.glob("event-*"))+
             sorted(root.glob("negative-*"))+sorted(root.glob("long-*")))
    result = [inspect(p) for p in paths if p.is_dir() and (p/"report.json").exists()]
    (root/"inspection.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    for item in result:
        print(json.dumps(item))
