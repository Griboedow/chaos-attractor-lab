"""Collect the verified cloud paths and compact numerical evidence."""
import json
from pathlib import Path
from search import CACHE, RESULTS, FAMILY, design
from phase_refine import BASES, local_design
from validate import write_json
import numpy as np


def collect():
    geometry=RESULTS/"geometry.json"
    geo={r["case"]:r for r in json.loads(geometry.read_text())["comparisons"]} if geometry.exists() else {}
    initials=[dict(np.load(p)) for p in sorted(CACHE.glob("batch-*/scout.npz"))]
    cycles=[json.loads(p.read_text()) for p in sorted(CACHE.glob("cycle-*/screen.json"))]
    hyper=[dict(np.load(p)) for p in sorted((CACHE/"hyper").glob("batch-*/scout.npz"))]
    nested=[json.loads(p.read_text()) for p in sorted((CACHE/"nested-phase").glob("screen-*.json"))]
    deep=[json.loads(p.read_text()) for p in sorted((CACHE/"nested-deep").glob("screen-*.json"))]
    ancestry={r["case"]:(data["base"],r["base_order"],r["multiplier"])
              for data in nested+deep for r in data["candidates"]}
    # Rebuild full batches only where the scout stored the passing subset.
    # Verify each saved passing row before including the reconstructed inputs.
    designs=[x["params"] for x in initials+hyper]
    for packet in cycles:
        batch=packet["batch"]
        params=design(batch,packet["count"])
        for hit in packet["hits"]:
            i=hit["case"]-1_000_000-batch*100_000
            if not np.array_equal(params[i],hit["params"]):
                raise ValueError("Mismatch in regenerated cyclic input")
        designs.append(params)
    for packet in nested+deep:
        base=packet["base"]
        group=next(i for i,(case,_) in enumerate(BASES) if case==base)
        row=json.loads((RESULTS/"checks"/f"{base}.json").read_text())
        params=local_design(np.array(row["params"]),group,packet["parameters"])
        for hit in packet["candidates"]:
            i=hit["case"]-20_000_000-group*100_000
            if not np.array_equal(params[i],hit["params"]):
                raise ValueError("Mismatch in regenerated local input")
        designs.append(params)
    all_inputs=np.vstack(designs)
    unique_inputs=len(np.unique(all_inputs,axis=0))
    summaries=[]
    rejected=[]
    for path in sorted((RESULTS/"checks").glob("*.json")):
        row=json.loads(path.read_text())
        if not row.get("cloud"):
            rejected.append(dict(case=row["case"],label=row["label"],reason=row["excluded"]))
            continue
        all_windows=min(row["reference_rates"])>0 and row["independent_tail"]["window_range"][0]>0
        phase_path=RESULTS/"phase"/path.name
        phase=json.loads(phase_path.read_text()) if phase_path.exists() else {}
        if not all_windows:
            rejected.append(dict(case=row["case"],label=row["label"],reason="Reference or independent leading windows included a non-positive value."))
            continue
        entry=dict(case=row["case"],label=row["label"],family=FAMILY,cloud=row["cloud"],key="points",
            cloud_relative=Path(row["cloud"]).relative_to(CACHE.parent).as_posix(),
            points=row["points"],sha256=row["sha256"],clearance_status="numerical candidate",
            mean_independent_spectrum=row["independent_spectrum"],
            positive_starts=sum(x["window_range"][0]>0 for x in row["basin"]),starts=len(row["basin"]),
            positive_nearby_parameters=sum(x["window_range"][0]>0 for x in row["local_parameters"]),
            nearby_parameters=len(row["local_parameters"]),check="checks/"+path.name)
        if phase.get("phase_separated"):
            entry.update(phase_groups=phase["order"],held_out_phase_accuracy=phase["held_out_accuracy"],
                         minimum_native_box_gap_fraction=phase["minimum_box_gap_fraction"])
            if "smallest_oriented_box_gap_fraction" in phase:
                entry["minimum_oriented_box_gap_fraction"]=phase["smallest_oriented_box_gap_fraction"]
                entry["separated_oriented_box_pairs"]=phase["pairs_with_positive_axis_gap"]
                entry["all_oriented_box_pairs"]=phase["finite_cloud_box_pairs"]
        long=RESULTS/"long-phase-controls"/path.name
        if long.exists():
            test=json.loads(long.read_text())
            entry["long_phase_controls"]="long-phase-controls/"+path.name
            entry["starts_with_all_1_5M_step_windows_positive"]=test["positive_all_long_windows"]
        f64=RESULTS/"phase-binary64"/path.name
        if f64.exists():
            test=json.loads(f64.read_text())
            entry["binary64_phase_controls"]="phase-binary64/"+path.name
            entry["combined_binary64_oriented_box_gap_fraction"]=test["combined_binary64"]["smallest_oriented_box_gap_fraction"]
        if row["case"] in geo:
            entry["nearest_archived_case"]=geo[row["case"]]["nearest_archived_case"]
            entry["nearest_archived_distance"]=geo[row["case"]]["archive_distance"]
        if row["case"] in ancestry:
            parent,order,multiple=ancestry[row["case"]]
            entry.update(parent_case=parent,parent_phase_order=order,phase_order_multiple=multiple)
        cert=RESULTS/"certificates"/f"{row['case']}-verification.json"
        if cert.exists():
            proof=json.loads(cert.read_text())
            entry["clearance_status"]="rigorous cyclic trapping; numerical chaotic motion"
            entry["certified_phase_count"]=proof["phase_regions"]
            entry["certificate"]="certificates/"+cert.name
            audit=cert.with_name(f"{row['case']}-directed.json")
            if audit.exists():
                entry["independent_certificate"]="certificates/"+audit.name
            uniform=cert.with_name(f"{row['case']}-parameter-neighborhood.json")
            if uniform.exists():
                entry["uniform_parameter_certificate"]="certificates/"+uniform.name
                entry["uniform_parameter_radius"]=2**-22
        native=RESULTS/"certificates"/f"{row['case']}-cycle-directed.json"
        if native.exists():
            audit=json.loads(native.read_text())
            entry["clearance_status"]="rigorous cyclic trapping; numerical chaotic motion"
            entry["certified_phase_count"]=audit["phase_regions"]
            entry["certificate"]="certificates/"+native.name
            second=native.with_name(f"{row['case']}-cycle-arb.json")
            if second.exists():
                entry["independent_certificate"]="certificates/"+second.name
        symbolic=RESULTS/"certificates"/f"{row['case']}-symbolic-directed.json"
        if symbolic.exists():
            audit=json.loads(symbolic.read_text())
            steps=audit["cover"]["map_steps"]
            entry["symbolic_certificate"]="certificates/"+symbolic.name
            entry["certified_entropy_lower_bound"]=f"log(2)/{steps} per original iterate"
            entry["lower_bound_value"]=float(np.log(2)/steps)
            entry["symbolic_set_scope"]="compact invariant subset of the global attractor"
            containment=RESULTS/"certificates"/f"{row['case']}-containment.json"
            if native.exists() and containment.exists():
                entry["symbolic_set_scope"]="compact invariant subset of the certified localized cyclic attracting set"
                entry["localization_check"]="certificates/"+containment.name
                entry["certified_least_periods"]=f"{steps} times every positive integer"
                for suffix in ("symbolic-arb","symbolic-arb-replay","periodic-points"):
                    name=f"{row['case']}-{suffix}.json"
                    extra=RESULTS/"certificates"/name
                    if extra.exists():
                        entry[suffix.replace("-","_")]="certificates/"+name
                uniform=RESULTS/"certificates"/f"{row['case']}-cycle-parameter-arb.json"
                if uniform.exists():
                    packet=json.loads(uniform.read_text())
                    entry["uniform_parameter_radius"]=2**-packet["parameter_radius_bits"]
                    entry["uniform_start_entry_step"]=packet["primary_start_entry"]["iterate"]
                    entry["uniform_trapping_certificate"]="certificates/"+uniform.name
                    entry["uniform_symbolic_certificate"]=f"certificates/{row['case']}-symbolic-parameter-arb.json"
                upper=RESULTS/"certificates"/f"{row['case']}-dimension-arb.json"
                if upper.exists():
                    entry["dimension_proof"]="dimension-proof.md"
                    entry["hausdorff_dimension_upper"]=json.loads(upper.read_text())["dimension_upper"]
                    lower=[json.loads((upper.parent/f"{row['case']}-dimension-phase-{kind}.json").read_text())["dimension_lower"]
                           for kind in ("arb","directed")]
                    entry["hausdorff_dimension_strict_lower_both_implementations"]=min(lower)
                    entry["topological_covering_dimension"]=1
                    roots=json.loads((upper.parent/f"{row['case']}-periodic-points.json").read_text())
                    entry["individually_isolated_saddle_orbits"]=len(roots["interval_certified_orbits"])
        precision=RESULTS/"precision"/path.name
        if precision.exists():
            probes=json.loads(precision.read_text())["probes"]
            entry["high_precision_leading_means"]={str(q["decimal_digits"]):q["mean"] for q in probes}
        if row["independent_spectrum"][1]>0:
            entry["independent_second_window_min"]=min(w[1] for w in row["independent_spectrum_windows"])
            two=RESULTS/"two-rate-controls"/path.name
            if two.exists():
                controls=json.loads(two.read_text())
                entry["starts_with_two_positive_rates_in_every_window"]=controls["all_second_windows_positive"]
        summaries.append(entry)
    initial_count=sum(len(r["case"]) for r in initials)
    initial_positive=sum(int((r["rates"].min(axis=1)>.012).sum()) for r in initials)
    initial_shaped=sum(int(np.isfinite(r["features"][:,0]).sum()) for r in initials)
    write_json(RESULTS/"manifest.json",dict(family=FAMILY,started_utc="2026-11-21",
        geometric_scout_designs=initial_count,geometric_scout_positive=initial_positive,
        geometric_scout_shape_evaluated=initial_shaped,phase_scout_designs=sum(c["count"] for c in cycles),
        phase_scout_filtered_positive=sum(c["screen_positive"] for c in cycles),
        phase_scout_hits=sum(len(c["hits"]) for c in cycles),
        two_exponent_scout_designs=sum(len(a["case"]) for a in hyper),
        two_exponent_scout_positive=sum(int((a["second_rates"].min(axis=1)>.009).sum()) for a in hyper),
        nested_phase_scout_designs=sum(a["parameters"] for a in nested),
        nested_phase_shortlist=sum(len(a["candidates"]) for a in nested),
        deeper_phase_scout_designs=sum(a["parameters"] for a in deep),
        deeper_phase_shortlist=sum(len(a["candidates"]) for a in deep),
        total_scout_parameter_rows=len(all_inputs),unique_parameter_tuples=unique_inputs,
        finalist_count=len(summaries),
        additional_basin_starts=4096*11,additional_basin_record="additional-basins.json",
        cloud_relative_root=str(CACHE.parent),cloud_relative_pattern="maps/clouds/<case>.npz",
        finalists=summaries,excluded=rejected))
    path=RESULTS/"README.md"
    text=path.read_text(encoding="utf-8")
    start,end="<!-- case-table-begin -->","<!-- case-table-end -->"
    if start in text and end in text:
        prefix=text.split(start,1)[0];suffix=text.split(end,1)[1]
        lines=["| Case and controls | Spatial form | Independent Lyapunov spectrum | Positive starts | Cyclic groups |",
               "| --- | --- | --- | --- | --- |"]
        for r in summaries:
            s=", ".join(f"{v:.4f}" for v in r["mean_independent_spectrum"])
            label=r.get("phase_groups","")
            if "certified_phase_count" in r:
                label=f"{r['certified_phase_count']} certified"
            lines.append(f"| [{r['case']}]({r['check']}) | {r['label']} | {s} | {r['positive_starts']}/{r['starts']} | {label} |")
        path.write_text(prefix+start+"\n\n"+"\n".join(lines)+"\n\n"+end+suffix,encoding="utf-8")
    print(len(summaries),"retained;",len(rejected),"excluded")


if __name__=="__main__":
    collect()
