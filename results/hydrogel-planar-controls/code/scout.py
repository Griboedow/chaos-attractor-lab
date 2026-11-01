"""Finite exploratory planar slice; does not test transverse 3D stability."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
import hydrogel_planar as model

ROOT = Path(__file__).resolve().parent
CASES = [(50,0.9,-33),(50,0.9,150),(50,0.9,300),(127,1.0,117),
         (127,1.0,200),(127,1.0,400),(200,1.0,160),(200,1.0,300)]

def run(item):
    i,(lam,eta,chi) = item
    folder = ROOT/"scout-01"/f"case-{i:02d}"
    folder.mkdir(parents=True,exist_ok=False)
    report = model.run_planar(32,model.parameters(slenderness=lam,eta=eta,chi=chi),
                              end=80,rtol=1e-6,output=folder/"trajectory.npz")
    report["model_sha256"] = hashlib.sha256(Path(model.__file__).read_bytes()).hexdigest()
    (folder/"report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    return i,{k:report[k] for k in ("parameters","success","elapsed_seconds","period_mean",
                                   "period_std","full_state_return_max_last10","angle_abs_max_late")}

if __name__ == "__main__":
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
        futures = [executor.submit(run,item) for item in enumerate(CASES)]
        for f in concurrent.futures.as_completed(futures):
            print(json.dumps(f.result()),flush=True)
