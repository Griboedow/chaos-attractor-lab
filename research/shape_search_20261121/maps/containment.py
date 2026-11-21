"""Confirm that the symbolic slab frame sits in the certified cyclic return box."""
import argparse
import json
from fractions import Fraction
from pathlib import Path


def check(inner,outer,dynamics):
    if inner["case"]!=outer["case"] or inner["params"]!=outer["params"]:
        raise ValueError("Different map parameters")
    if inner["center"]!=outer["center"] or inner["basis"]!=outer["basis"]:
        raise ValueError("This exact comparison requires the same affine frame")
    if dynamics["cover"]["map_steps"]!=outer["case_order"]:
        raise ValueError("The symbolic and cyclic returns use different iterates")
    q=Fraction.from_float
    low=[q(i)-q(o) for i,o in zip(inner["lo"][0],outer["lo"][0])]
    high=[q(o)-q(i) for i,o in zip(inner["hi"][0],outer["hi"][0])]
    if min(low+high)<0:
        raise ValueError("Symbolic box is not contained in the trapping box")
    cover=dynamics["cover"]
    if not (-1<=cover["faces"][0]<cover["faces"][1]<cover["faces"][2]<cover["faces"][3]<=1):
        raise ValueError("The symbolic slabs leave their frame")
    return dict(case=outer["case"],return_order=outer["case_order"],
                symbolic_box_contained=True,lower_bound_differences_exact=list(map(str,low)),
                upper_bound_differences_exact=list(map(str,high)),
                consequence="The full symbolic subsystem lies in the localized cyclic attracting set.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("inner","outer","symbolic"):
        p.add_argument(k,type=Path)
    p.add_argument("--output",type=Path,required=True);a=p.parse_args()
    result=check(*(json.loads(getattr(a,k).read_text()) for k in ("inner","outer","symbolic")))
    a.output.write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)
