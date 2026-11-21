"""Count parameter vectors separately from replays, continuations and output clouds."""
import argparse
import csv
import importlib.util
import math
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np
from index_campaign import read,write,RESULTS

HERE=Path(__file__).resolve().parent


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def audit(cache):
    braid=module(HERE/"braids/search.py","fourier_designs")
    receipt=read(RESULTS/"braids/certificates.json")
    rng=np.random.default_rng(receipt["seed"])
    braid_keys=set()
    for i in range(receipt["draws"]):
        row=braid.draw_case(rng,i+1)
        coeff=row["harmonics"]
        # Treat scalar multiples as potentially equal before the data-dependent
        # scale is known. This count is conservative for the actual screened maps.
        common=math.gcd(*(int(c) for k,x,y in coeff for c in (x,y)))
        waves=tuple(sorted((k,x//common,y//common) for k,x,y in coeff))
        braid_keys.add((row["degree"],row["fiber_turns"],waves))
    pinch=module(HERE/"pinched/search.py","pinched_designs")
    pinched={ (r["kappa"],tuple(sorted(map(tuple,r["waves"])))) for r in pinch.cases(1536)}
    delay_keys=set()
    delay_rows=0
    for group in ["scout","recruits"]:
        with (cache/"delay"/group/"screen.csv").open(newline="",encoding="utf-8") as f:
            for row in csv.DictReader(f):
                values=[float(row[k]) for k in ("r","a","b","c")]
                if "survival" in row:
                    values.append(float(row["survival"]))
                delay_keys.add(tuple(values));delay_rows+=1
    # The local screen writes only survivors. Replay its seeded parameter draws
    # with its numerical evaluator replaced by an observer.
    sys.path.insert(0,str(HERE/"delay"))
    local=module(HERE/"delay/neighborhood.py","delay_designs")
    def observe(par,*args):
        nonlocal delay_rows
        delay_rows+=1;delay_keys.add(tuple(map(float,par)))
        raise ValueError("Observer has recorded this parameter vector")
    with patch.object(local,"native_spectrum",observe):
        local.run(cache/"delay/recruits/screen.json",cache/"delay/parameter-audit",32000)
    twist=read(RESULTS/"maps/manifest.json")
    proj=read(RESULTS/"projection/distinct-configurations.json")
    flows=read(RESULTS/"flows/index-summary.json")
    assert delay_rows==104000
    counts={"projection":proj["distinct_digests"],"twist":twist["unique_parameter_tuples"],"flows":flows["unique_screened_equations"],
            "delay":len(delay_keys),"fourier_braids_lower_bound":len(braid_keys),"pinched":len(pinched)}
    write(RESULTS/"design-counts.json",{"designs_by_track":counts,"total_lower_bound":sum(counts.values()),
        "delay_rows_attempted":delay_rows,"braid_rows_attempted":receipt["draws"],"pinched_rows_attempted":1536,
        "projection_receipt":"projection/distinct-configurations.json",
        "twist_receipt":"maps/manifest.json","flow_receipt":"flows/index-summary.json",
        "scope":"Broad and targeted parameter screens. Continuations, replays and nonautonomous phase observations are excluded."})
    print(counts,sum(counts.values()),flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache",required=True,type=Path)
    audit(ap.parse_args().cache)
