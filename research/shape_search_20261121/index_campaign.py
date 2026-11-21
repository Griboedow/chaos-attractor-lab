"""Build the campaign indexes from completed, independently replayable experiments."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parents[2] / "results/shape_search_20261121"


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path, data):
    path.parent.mkdir(parents=True,exist_ok=True)
    content = (json.dumps(data, indent=2, allow_nan=False) + "\n").encode()
    path.write_bytes(gzip.compress(content, mtime=0) if path.suffix == ".gz" else content)


def receipt(row, folder, cache):
    path = folder / row["cloud"]
    with (cache / path).open("rb") as f:
        sha = hashlib.file_digest(f,"sha256").hexdigest()
    return {"case":row["case"], "cloud_relative":path.as_posix(), "file_sha256":sha}


def index(cache):
    inputs = {}
    rows = []
    for group in ("scout-long","recruits-long","local-long","entropy-box"):
        folder = Path("delay") / group
        for run in read(cache / folder / "manifest.json"):
            item = receipt(run,folder,cache)
            item.update({k:v for k,v in run.items() if k!="cloud"})
            item["positive_primary"] = min(v[0] for v in item["native_spectra"])>0
            rows.append(item)
    write(RESULTS/"delay/manifest.json",rows)
    for name in ["delay","braids","pinched"]:
        if name == "pinched":
            folder = Path("pinched/clouds-v2")
            samples = read(cache / folder / "manifest.json")
            rows = [{**r,**receipt(r,folder,cache)} for r in samples]
            write(RESULTS / "pinched/manifest.json",{"count_screened":1536,"selected":rows})
        if name == "braids":
            folder = Path("braids/clouds")
            data = read(cache / folder / "manifest.json")
            rows = [{**receipt(r,folder,cache),"degree": r["degree"],
                     "core_determinant":r["determinant"],
                     "base_exponents":[v["base_exponent"] for v in r["runs"]]} for r in data["rows"]]
            write(RESULTS / "braids/manifest.json",
                  {"steps_per_initial":data["steps_per_initial"],"burn":data["burn"],"rows":rows})
        inputs.update({f"{name}/{r['case']}":r["cloud_relative"] for r in rows})
    projections = read(RESULTS / "projection/manifest.json")
    for file in sorted((RESULTS/"projection/three").glob("case-*.json")):
        run=read(file)
        if run.get("three_positive_passes"):
            cloud=Path(run["cloud"]).relative_to(cache).as_posix()
            inputs[f"projection/{run['case']}"]=cloud
    for file in sorted((RESULTS/"projection/three/companions").glob("companion-*.json")):
        run=read(file)
        inputs[f"companions/{run['case']}"]=Path(run["cloud"]).relative_to(cache).as_posix()
    for run in read(RESULTS/"delay/sink-clouds.json"):
        inputs[f"delaycycles/{run['period']}"]=run["cloud_relative"]
    maps = read(RESULTS / "maps/manifest.json")["finalists"]
    flow_data = read(RESULTS / "flows/manifest.json")
    for name,rows in [("projection",projections),("maps",maps),("flows",flow_data["cases"]),
                      ("cycles",flow_data["stable_cycles"]),("links",flow_data["linked_cycles"]),
                      ("skeletons",flow_data["knotted_periodic_skeletons"])]:
        inputs.update({f"{name}/{r['case']}":r["cloud_relative"] for r in rows})
    with (RESULTS/"flows/periodic-catalog.csv").open(newline="",encoding="utf-8") as f:
        for cycle in csv.DictReader(f):
            inputs.setdefault(f"cycles/{cycle['case']}",cycle["cloud_relative"])
    for orbit in read(RESULTS/"flows/nested-skeletons.json")["cores"]:
        inputs[f"skeletons/{orbit['case']}-order-{orbit['period_order']}"]=orbit["cloud_relative"]
    companion = read(RESULTS/"flows/final-link-crosscheck.json")["symmetric_companion"]
    inputs["skeletons/180844-order-2-symmetric"]=companion["cloud_relative"]
    pairs = read(RESULTS/"presentation.json")
    selected = {}
    for pair in pairs:
        for item in pair["cases"]:
            track=item.get("track",pair["track"])
            selected[f"{track}/{item['case']}"] = item["cloud"]
    write(RESULTS/"comparison-inputs.json",inputs)
    write(RESULTS/"selected-comparison-inputs.json",selected)
    refs = read(cache/"extended-reference.json")
    write(RESULTS/"reference-receipts.json",{"count":len(refs["rows"]),"references":[
        {k:v for k,v in r.items() if k!="descriptors"} for r in refs["rows"]]})
    write(RESULTS/"delay/period-39-basins.json", read(cache/"delay/period39-basins/summary.json"))
    return len(inputs),len(selected)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache",required=True,type=Path)
    print(index(ap.parse_args().cache))
