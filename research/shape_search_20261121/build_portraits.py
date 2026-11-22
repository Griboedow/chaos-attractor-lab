"""Register selected run inputs and build their canonical portraits and rotations."""
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import animate_portraits
import render_portraits
from data_cache import cache_root, cached_path, digest_file

RESULT = ROOT/"results/shape_search_20261121"
OUT = ROOT/"docs/assets/portraits"


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp = path.with_name(path.name+".tmp")
    temp.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n",encoding="utf-8",newline="\n")
    temp.replace(path)


def register_input(raw, virtual, registry, group, command):
    checksum = digest_file(raw)
    old = next((r for r in registry["files"] if r["path"]==virtual),None)
    if old:
        assert old["sha256"]==checksum, (virtual,"Existing input changed")
    else:
        with np.load(raw) as f:
            arrays = [{"key":k,"shape":list(f[k].shape),"dtype":str(f[k].dtype),
                       "bytes":f[k].nbytes,"retained":False} for k in f.files]
        old = {"path":virtual,"bytes":raw.stat().st_size,"git_blob":None,"sha256":checksum,
               "group":group,"compact":None,"arrays":arrays,"recovery":"recompute",
               "recompute_command":command}
        registry["files"].append(old)
        registry["removed_bytes"] += old["bytes"]
    dest = cached_path(old)
    if dest.exists():
        assert digest_file(dest) == checksum
    else:
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(raw,dest)
        assert digest_file(dest) == checksum


def make_row(study,runs,registry):
    ident = study["id"]
    default_projection = render_portraits.camera(*study["view"]).tolist()
    panels,animation = [],[]
    for case in study["cases"]:
        number = str(case["case"])
        virtual = f"results/shape_search_20261121/clouds/{study['track']}-{number}.npz"
        raw = runs/case["cloud"]
        register_input(raw,virtual,registry,"shape-search-20261121",case["reproduce"])
        # The atlas renderer supports archived count grids. Define the physical
        # camera here and use the same projection for both renderers.
        with np.load(raw) as f:
            z = np.asarray(f["points"],dtype=np.float64)
        projection=default_projection
        if study.get("pca_camera"):
            centered=z-z.mean(0)
            _,axes=np.linalg.eigh(np.cov(centered.T))
            axes=axes[:,::-1]
            skew=np.mean((centered@axes)**3,axis=0)
            assert np.all(np.abs(skew)>1e-12)
            axes=axes*np.sign(skew)
            s=np.sqrt(3)/2
            frame=np.array([[s,0,.5],[-.25,s,s/2]])
            projection=(frame@axes.T).tolist()
        center = (z.min(0)+z.max(0))/2
        xy = (z-center)@np.asarray(projection).T
        radius = float(np.abs(xy).max())*1.06
        counts = np.histogram2d(xy[:,0],xy[:,1],bins=760,range=[[-radius,radius]]*2)[0].T[::-1].astype(np.uint32)
        assert counts.sum()==len(z)
        counts_path = runs/"portrait-counts"/(ident+"-"+number+".npz")
        counts_path.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(counts_path,counts=counts,projection=projection,center=center,half_width=radius)
        counts_name = f"results/shape_search_20261121/clouds/{counts_path.stem}-{digest_file(counts_path)[:12]}.npz"
        register_input(counts_path,counts_name,registry,"shape-search-20261121",
                       "Regenerate the raw run, then run python research/shape_search_20261121/build_portraits.py.")
        panels.append({"label":number,"source":counts_name,"key":"counts","mode":"histogram"})
        meta_path = RESULT/"cameras"/(ident+"-"+number+".json")
        write(meta_path,{"center":center.tolist(),"basis":np.eye(3).tolist(),"projection":projection,
                         "source":virtual,"key":"points","points":len(z),"grid_size":760,
                         "half_width":radius,"source_sha256":digest_file(raw)})
        if study["featured"]:
            n=len(z) if study.get("all_components") else min(len(z),60000)
            animation.append({"source":virtual,"key":"points","last_count":n,
                              "metadata_source":meta_path.relative_to(ROOT).as_posix()})
    row = {"id":ident,"title":study["title"],"group":study["group"],
           "report":"docs/shape-search-20261121.md","source_portrait":f"docs/assets/portraits/{ident}.png",
           "description":study["description"],"coordinates":study["coordinates"],
           "featured":study["featured"],"projection":default_projection,"panels":panels}
    if animation:
        row["animation_panels"]=animation
    return row


def update_readme(catalog):
    path = ROOT / "README.md"
    main = path.read_text(encoding="utf-8")
    studies = [row for row in catalog["rows"] if row["id"].startswith("shapes-20261121-")]
    for featured, position in ((True, "Further portraits"), (False, "Research index")):
        name = "selected" if featured else "additional"
        start = f"<!-- shape-search-20261121-{name} -->"
        end = f"<!-- /shape-search-20261121-{name} -->"
        heading = "## Shape search of 2026-11-21" if featured else "### Spatial maps and flows"
        title_level = "###" if featured else "####"
        entries = []
        for row in studies:
            if bool(row.get("featured")) != featured:
                continue
            title = row["title"]
            stem = "docs/assets/portraits/" + row["id"]
            entries.extend([f"{title_level} {title}", row["description"]])
            if featured:
                entries.append(f"![{title} in rotation]({stem}.webp)")
            entries.extend([
                f"![{title}]({stem}.png)",
                row["coordinates"] + " [Equations, controls and prior work](docs/shape-search-20261121.md)",
            ])

        block = start + "\n" + "\n\n".join([heading, *entries]) + "\n" + end if entries else ""
        if start in main:
            first = main.index(start)
            last = main.index(end, first) + len(end)
            main = main[:first] + block + main[last:]
        elif block:
            anchor = f"## {position}\n"
            assert main.count(anchor) == 1
            main = main.replace(anchor, block + "\n\n" + anchor)
    path.write_text(main, encoding="utf-8", newline="\n")
    subprocess.run([sys.executable, str(ROOT / "tools/render_dates.py"), "sync"], check=True)
    subprocess.run([sys.executable, str(ROOT / "tools/portrait_equations.py"), "sync"], check=True)


def build(specs,runs,refresh=False):
    registry_path = ROOT/"records/generated-data.json"
    catalog_path = ROOT/"records/portrait-catalog.json"
    registry,catalog = read(registry_path),read(catalog_path)
    rows = []
    for study in specs:
        row = make_row(study,runs,registry)
        previous = next((s for s in catalog["rows"] if s["id"]==row["id"]),None)
        if previous is None:
            catalog["rows"].append(row)
        elif refresh:
            catalog["rows"][catalog["rows"].index(previous)] = row
        else:
            text_fields={"title","group","description","coordinates","report"}
            assert {k:v for k,v in row.items() if k not in text_fields}=={k:v for k,v in previous.items() if k not in text_fields},row["id"]
            previous.update({k:row[k] for k in text_fields})
        rows.append(row)
    write(registry_path,registry)
    write(catalog_path,catalog)
    # The existing manifests are cumulative. Preserve each earlier image receipt.
    render_manifest = read(OUT/"render.json")
    anim_manifest = read(OUT/"animation.json")
    if refresh:
        ids = {row["id"] for row in rows}
        render_manifest["rows"] = [r for r in render_manifest["rows"] if r["id"] not in ids]
        anim_manifest["rows"] = [r for r in anim_manifest["rows"] if r["id"] not in ids]
    for row in rows:
        if not any(r["id"] == row["id"] for r in render_manifest["rows"]):
            receipt=render_portraits.render(row,OUT)
            subprocess.run([sys.executable,str(ROOT/"tools/render_dates.py"),
                            "record",receipt["file"]],check=True)
            receipt.update(renderer_sha256=digest_file(ROOT/"tools/render_portraits.py"),
                           data_resolver_sha256=digest_file(ROOT/"tools/data_cache.py"))
            render_manifest["rows"].append(receipt)
            render_manifest["catalog_sha256"]=digest_file(catalog_path)
            write(OUT/"render.json",render_manifest)
        if row["featured"] and not any(r["id"] == row["id"] for r in anim_manifest["rows"]):
            receipt=animate_portraits.render(row,OUT)
            receipt.update(renderer_sha256=digest_file(ROOT/"tools/animate_portraits.py"),
                           static_renderer_sha256=digest_file(ROOT/"tools/render_portraits.py"),
                           data_resolver_sha256=digest_file(ROOT/"tools/data_cache.py"))
            anim_manifest["rows"].append(receipt)
            anim_manifest["catalog_sha256"]=digest_file(catalog_path)
            write(OUT/"animation.json",anim_manifest)
    render_manifest["catalog_sha256"]=anim_manifest["catalog_sha256"]=digest_file(catalog_path)
    order = {row["id"]:n for n,row in enumerate(catalog["rows"])}
    for book in (render_manifest,anim_manifest):
        book["rows"].sort(key=lambda r:order[r["id"]])
    write(OUT/"render.json",render_manifest)
    write(OUT/"animation.json",anim_manifest)
    update_readme(catalog)


if __name__ == "__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs",type=Path,default=cache_root()/"shape-search-20261121")
    ap.add_argument("--selection",type=Path,default=RESULT/"presentation.json")
    ap.add_argument("--refresh",action="store_true",help="Render the selected rows again after a camera change.")
    args=ap.parse_args()
    build(read(args.selection),args.runs,args.refresh)
