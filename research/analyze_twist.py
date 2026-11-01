"""Summarize complete cohorts; retain regular outcomes and sampling limitations."""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from twist_map import step_tangent


def main():
    base=Path('results/twist-geometry')
    long=np.load(base/'long-runs.npz')
    basins=np.load(base/'basins.npz')
    result={'candidates':[]}
    for i,ident in enumerate(long['ids']):
        spectrum=long['windows'][i].mean(0)
        ordered=np.sort(spectrum)[::-1]
        j=int((np.cumsum(ordered)>0).sum())
        dimension=j+ordered[:j].sum()/abs(ordered[j])
        subset=np.where(basins['ids']==ident)[0]
        positive=basins['windows'][subset,:,0].min(1)>0
        tree=cKDTree(long['tail'][i])
        diameter=np.linalg.norm(np.ptp(long['tail'][i],axis=0))
        distances=[];regular=[]
        for lane,is_positive in zip(subset,positive):
            z=basins['tail'][lane]
            if is_positive:
                distances.append(float(np.quantile(tree.query(z,workers=1)[0],.99)/diameter))
            else:
                record={'lane':int(lane),'period':None}
                for period in range(1,2049):
                    residual=float(np.max(np.linalg.norm(z[period:]-z[:-period],axis=1)))
                    if residual<1e-11:
                        state=z[-1].copy();jac=np.eye(3)
                        for _ in range(period):
                            state,jac=step_tangent(state,jac,basins['parameters'][lane])
                        record.update(period=period,max_tail_return_error=residual,
                                      monodromy_spectral_radius=float(np.max(np.abs(np.linalg.eigvals(jac)))))
                        break
                regular.append(record)
        row={'id':int(ident),'spectrum':spectrum.tolist(),'kaplan_yorke_estimate':float(dimension),
             'leading_window_range':[float(long['windows'][i,:,0].min()),float(long['windows'][i,:,0].max())],
             'max_volume_identity_error':float(np.max(np.abs(long['windows'][i].sum(1)-3*np.log(long['parameters'][i,0])))),
             'max_observed_radius':float(long['radii'][i]),'basin_positive_all_windows':int(positive.sum()),
             'basin_tested':len(subset),'positive_basin_tail_distance_99pct_max_over_bbox_diagonal':max(distances),
             'regular_basin_outcomes':regular}
        result['candidates'].append(row)
    if (base/'neighbors.npz').exists():
        d=np.load(base/'neighbors.npz');result['neighbors']=[]
        for ident in np.unique(d['ids']):
            select=d['ids']==ident;w=d['windows'][select,:,0]
            result['neighbors'].append({'id':int(ident),'tested':int(select.sum()),
                                        'positive_all_windows':int((w.min(1)>0).sum()),
                                        'mean_exponent_range':[float(w.mean(1).min()),float(w.mean(1).max())],
                                        'rounded_to_3_decimals_spectrum':d['windows'][select][-1].mean(0).tolist()})
    (base/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
