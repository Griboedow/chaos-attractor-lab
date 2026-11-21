"""Write the measured search totals and nearest-reference controls into the report."""
import gzip
import json
from pathlib import Path
import re

import numpy as np
from index_campaign import read,write,RESULTS


def cached(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def update():
    cloud=read(RESULTS/"cloud-totals.json")
    design=read(RESULTS/"design-counts.json")
    all_results=cached(RESULTS/"comparison.json.gz")
    selected=cached(RESULTS/"selected-comparison.json.gz")
    assert len(all_results)==cloud["clouds"]
    assert set(all_results)==set(read(RESULTS/"comparison-inputs.json"))
    assert set(selected)==set(read(RESULTS/"selected-comparison-inputs.json"))
    shapes=cloud["by_group"]
    selection={}
    for k,v in selected.items():
        row={}
        for mode,check in v["compared"].items():
            row[mode]={"nearest_reference":check["nearest_reference"][0][1],
                       "fingerprint_distance":check["nearest_reference"][0][0],
                       "split_distance":check["split_distance"],
                       "registration_residual":check["local_icp_to_nearest"],
                       "split_registration_residual":check["local_icp_between_halves"]}
        selection[k]=row
    summary={}
    for mode in ["similarity","affine"]:
        refs=set()
        for k,r in all_results.items():
            if mode in r["compared"]:
                refs.add(r["compared"][mode]["reference_count"])
        ratios=[d[mode]["registration_residual"]/d[mode]["split_registration_residual"]
                for d in selection.values() if mode in d and d[mode]["split_registration_residual"]>0]
        summary[mode]={"reference_counts":sorted(refs),"aligned_selected":len(ratios),
                       "alignment_ratio_quartiles":np.quantile(ratios,[.25,.5,.75]).tolist(),
                       "ratio_exceeds_one":sum(v>1 for v in ratios)}
    write(RESULTS/"comparison-summary.json",{"clouds":len(all_results),"by_mode":summary,"selected":selection})
    lines=["## Search totals","",
           f"The campaign screened at least **{design['total_lower_bound']:,} distinct parameter settings** across six search tracks. The archive contains **{cloud['clouds']:,} checked three-coordinate samples** with {cloud['points']:,} stored points. It includes geometry controls, periodic skeletons and multiple coexisting sets at the same parameters. All 30 gallery input files were regenerated from the recorded parameters and matched their archived SHA-256 checksums byte for byte. [Design audit](../results/shape_search_20261121/design-counts.json) · [Cloud counts](../results/shape_search_20261121/cloud-totals.json) · [Gallery replay](../results/shape_search_20261121/gallery-replay.json) · [Independent array checks](../results/shape_search_20261121/replay-audit.json)",
           "","| Search family | Screened settings | Archived samples |","| --- | ---: | ---: |"]
    groups=[("Toroidal projections","projection",["projection","companions"]),("Three-axis twists","twist",["maps"]),
            ("Continuous flows","flows",["flows","cycles","links","skeletons"]),
            ("Fourier solid-torus maps","fourier_braids_lower_bound",["braids"]),
            ("Delayed population maps","delay",["delay","delaycycles"]),
            ("Pinched graphs","pinched",["pinched"])]
    for label,key,subgroups in groups:
        n=design["designs_by_track"][key]
        val=sum(shapes[g]["clouds"] for g in subgroups)
        lines.append(f"| {label} | {n:,} | {val:,} |")
    lines.extend(["","The Fourier count merges proportional coefficient lists before normalization and gives a lower bound. Local sweeps used solely to locate bifurcations and replicate trajectories are excluded from the design total."])
    doc=RESULTS.parents[1]/"docs/shape-search-20261121.md"
    s=doc.read_text(encoding="utf-8")
    start,end="<!-- campaign-totals -->","<!-- /campaign-totals -->"
    section=start+"\n"+"\n".join(lines)+"\n"+end+"\n\n"
    s=re.sub(re.escape(start)+".*?"+re.escape(end)+"\n\n","",s,flags=re.S)
    anchor="## Follow-up search" if "## Follow-up search" in s else "## A knotted chaotic band"
    s=s.replace(anchor,section+anchor,1)
    doc.write_text(s,encoding="utf-8",newline="\n")
    return summary


if __name__=="__main__":
    print(update())
