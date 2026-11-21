"""Analyze consecutive pairs with the same frozen halfspace label in binary64 traces."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from capture import initial_state


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf8")


def median_range(values):
    if not len(values):
        return []
    return np.quantile(values, [0, .1, .5, .9, 1]).tolist()


def fit(z):
    normal = z[1::2].mean(0) - z[::2].mean(0)
    normal /= np.linalg.norm(normal)
    cut = float(((z[::2] @ normal).max() + (z[1::2] @ normal).min()) / 2)
    return normal, cut


def make_planes(source, saved):
    n, c = np.asarray(saved["normal"]), saved["cut"]
    center, rms = source.mean(0), saved["source_rms"]
    pivot = center + (c - center @ n) * n
    u = np.cross(n, [1, 0, 0])
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    planes = [("frozen", n, c), ("first-half fit", *fit(source[:len(source)//2]))]
    for offset in [-.25, -.1, .1, .25]:
        planes.append((f"offset {offset:+.2f} RMS", n, c + offset * rms))
    for label, axis in (("u", u), ("v", v)):
        for degrees in (-7, 7):
            theta = np.deg2rad(degrees)
            normal = n * np.cos(theta) + axis * np.sin(theta)
            planes.append((f"tilt {label} {degrees:+d} deg", normal, pivot @ normal))
    return planes


def event_times(z, normal, cut):
    ids = []
    for start in range(0, len(z) - 1, 1_000_000):
        block = z[start:min(start + 1_000_001, len(z))]
        side = block @ normal > cut
        ids.append(start + np.flatnonzero(side[1:] == side[:-1]))
    return np.concatenate(ids)


def runs(events, max_gap):
    pieces = np.split(events, np.flatnonzero(np.diff(events) > max_gap) + 1) if len(events) else []
    return [(int(p[0]), int(p[-1]), len(p)) for p in pieces]


def summarize_runs(events, gap):
    clusters = runs(events, gap)
    return dict(clusters=len(clusters), odd=sum(n % 2 for _, _, n in clusters),
                maximum_edges=max((n for _, _, n in clusters), default=0),
                maximum_span=max((b - a + 1 for a, b, _ in clusters), default=0))


def geometry(x, p):
    sites = np.array(p["sites"])
    weights = -p["inv_tau"] * np.sum((x[:, None, :] - sites[None, :, :]) ** 2, axis=2)
    weights = np.exp(weights - weights.max(1)[:, None])
    weights /= weights.sum(1)[:, None]
    h = (2 * (weights @ sites) - x) / np.array(p["scales"])
    h[:, 2] -= p["warp"] * (h[:, 0] ** 2 - h[:, 1] ** 2)
    rho = np.linalg.norm(h[:, :2], axis=1)
    radius = np.hypot(rho - 1, h[:, 2])
    return np.column_stack((rho, radius, weights.max(1)))


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def main(args):
    meta = json.loads(args.reference.read_text())
    doc = json.loads(args.record.read_text())
    with np.load(args.cloud) as f:
        source = f["points"].astype("float64")
    planes = make_planes(source, meta)
    normals, cuts = {}, {}
    for label, normal, cut in planes:
        normals[label], cuts[label] = normal.tolist(), float(cut)
    rms = meta["source_rms"]
    normal, cut = planes[0][1:]
    ref_tree = cKDTree(source)
    args.output.mkdir(parents=True, exist_ok=True)
    for path in sorted(args.cache.glob("trace-*.npy")):
        if ".partial." in path.name:
            continue
        cohort = int(path.name.split("-")[1])
        summary = args.output / f"cohort-{cohort:02d}.json"
        if summary.exists():
            continue
        z = np.load(path, mmap_mode="r")
        times = event_times(z, normal, cut)
        if cohort == 0:
            assert len(z) == 50000000 and len(times) == meta["tests"][0]["same_halfspace_successive_pairs"]
            assert np.array_equal(z[4750000:5000000].astype("f4"), source.astype("f4"))
        assert not len(times) or times.min() >= 16 and times.max() + 17 < len(z)
        offsets = np.arange(-16, 18)
        windows = z[times[:, None] + offsets]
        src = z[times]
        radii = geometry(src, doc["parameters"]) if len(src) else np.empty((0, 3))
        labels = src @ normal > cut
        ctrl_ids = np.random.default_rng(855476 + cohort).integers(0, len(z), 10000)
        ctrl = z[ctrl_ids]
        ctrl_geo = geometry(ctrl, doc["parameters"])
        # Physical distances to saved image points measure departures from both lobes.
        ref_dist = ref_tree.query(windows.reshape(-1, 3))[0].reshape(len(windows), len(offsets)) / rms
        control_dist = ref_tree.query(ctrl)[0] / rms
        variants = {}
        clusters = runs(times, 64)
        for label, n, delta in planes:
            other = times if label == "frozen" else event_times(z, n, delta)
            shared = int(np.intersect1d(times, other).size)
            within = np.zeros(len(times), dtype=bool)
            if len(other) and len(times):
                i = np.searchsorted(other, times)
                within = np.minimum(np.abs(times-other[i.clip(0,len(other)-1)]), np.abs(times-other[(i-1).clip(0,len(other)-1)])) <= 4
            variants[label] = dict(count=len(other), shared=shared, original_with_match_within_4=int(within.sum()),
                                    clusters_64=summarize_runs(other, 64)["clusters"])
        blocks = np.bincount(times // 1000000, minlength=(len(z) + 999999) // 1000000)
        down = np.abs((windows @ normal - cut) / rms)
        outfile = args.cache / f"event-windows-{cohort:02d}.npz"
        np.savez_compressed(outfile, left_index=times, offsets=offsets, xyz=windows, pre_geometry=radii,
                            cloud_distance_rms=ref_dist)
        distance_q = [median_range(a) for a in radii.T]
        control_q = [median_range(a) for a in ctrl_geo.T]
        write(summary, dict(cohort=cohort, observed=len(z), compared=len(z)-1, events=len(times),
              initial=initial_state(doc,cohort).tolist(), first=int(times.min()) if len(times) else None,
              last=int(times.max()) if len(times) else None, million_step_counts=blocks.tolist(),
              side_positive=int(labels.sum()), side_negative=int((~labels).sum()),
              inter_event_steps=median_range(np.diff(times)), groupings={str(k):summarize_runs(times,k) for k in [1,16,64,256,4096]},
              geometry_columns=["axis_radius", "core_distance", "peak_posterior"], geometry_event_quantiles=distance_q,
              geometry_control_quantiles=control_q, quantiles=[0,.1,.5,.9,1],
              normal=normal.tolist(), cut=cut, rms=rms, variants=variants,
              nearest_cloud_control=median_range(control_dist), nearest_cloud_event=median_range(ref_dist[:,16]),
              nearest_cloud_at_plus_16=median_range(ref_dist[:,-1]),
              min_endpoint_plane_margin=median_range(np.min(down[:,16:18],axis=1)),
              window_offsets=[-16,17], trace=path.name, windows=outfile.name, windows_sha256=sha(outfile),
              events_relative_to_first_recorded_state="zero-based index of left endpoint; first stored state follows discarded burn and one map iteration"))
        print(cohort,len(times),len(clusters),flush=True)
    write(args.output/"planes.json", dict(normals=normals, cuts=cuts, cloud_sha256=sha(args.cloud), rms=rms,
          holdout_method="fit first half, inspect second half",
          heldout_first_half_classifier=dict(pairs=len(source)//2-1, same_halfspace=len(event_times(source[len(source)//2:],*planes[1][1:])),
          negative_margin_rms=float((planes[1][2]-(source[len(source)//2::2]@planes[1][1]).max())/rms),
          positive_margin_rms=float(((source[len(source)//2+1::2]@planes[1][1]).min()-planes[1][2])/rms))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for name in ("record","reference","cloud","cache","output"):
        parser.add_argument("--"+name,required=True,type=Path)
    main(parser.parse_args())
