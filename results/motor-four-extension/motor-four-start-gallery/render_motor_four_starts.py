"""Compare the eight observed symmetry classes without merging their trajectories."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from audit_motor_permutation_geometry import full
from render_motor_four import H
from render_motor_four_ribbon import style


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf8')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    report=json.loads((args.source/'comparisons.json').read_text())
    assert report['completed']==report['declared']==32
    classes={};runs=[];counts=Counter()
    for record in report['runs']:
        row=record['distinguished_by_smallest_mean_position']
        sign=int(np.sign(record['circulation_of_sorted_other_rows']))
        label=(row,sign);counts[label]+=1;classes.setdefault(label,record)
        best=record['best_by_median']
        near=[r['permutation'] for r in record['comparisons'] if r['median']<=3*best['median'] and r['p95']<=3*best['p95']]
        runs.append(dict(id=record['id'],distinguished=row,circulation_sign=sign,best=best,near_permutations=near))
    assert len(classes)==8
    records=[classes[label] for label in sorted(classes)]
    arrays=[];positions=[];aligned=[];hashes=[]
    for record in records:
        path=args.source/record['file'];digest=hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest==record['sha256']
        data=np.load(path);ids=np.flatnonzero(data['time']>=data['time'][-1]-600)[::2]
        x=data['states'][ids];position=full(x,4)[:,8:]
        arrays.append(x);positions.append(position@H.T)
        aligned.append(position[:,record['best_by_median']['permutation']]@H.T)
        hashes.append(dict(id=record['id'],source=str(path),sha256=digest,display_samples=len(ids),
            display_start=float(data['time'][ids[0]]),display_end=float(data['time'][ids[-1]])))
    np.savez_compressed(out/'display-inputs.npz',states=np.array(arrays),ids=[r['id'] for r in records],
        permutations=[r['best_by_median']['permutation'] for r in records],contrast=H)
    bg='#071320';fg='#e9f0f5';muted='#99afc0'
    colors=['#4bd4cd','#f4b879','#a2a8ff','#ee82a4','#72b4eb','#cadb80','#cf91e8','#f0cd8d']
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg,'axes.labelcolor':muted,
        'xtick.color':muted,'ytick.color':muted})
    allpoints=np.concatenate(positions)
    fig=plt.figure(figsize=(16,10.5),facecolor=bg)
    for k,(record,points) in enumerate(zip(records,positions)):
        ax=fig.add_subplot(2,4,k+1,projection='3d',facecolor=bg)
        ax.plot(*points.T,lw=.32,alpha=.67,color=colors[k])
        style(ax,allpoints,27,-45,zoom=1.03)
        key=(record['distinguished_by_smallest_mean_position'],int(np.sign(record['circulation_of_sorted_other_rows'])))
        ax.set_title(f'Row {key[0]} distinguished | orientation {key[1]:+d}\nStart {record["id"]:02d} | {counts[key]} of 32 starts',fontsize=11,pad=0)
    fig.suptitle('ONE FAMILY, EIGHT OBSERVED ORIENTATIONS',x=.04,ha='left',y=.98,fontsize=23,weight='bold')
    fig.text(.04,.915,'Each panel is a separately integrated fresh start, with the same physical coordinates, scale and camera.\nOne representative per observed class; final 600 of 3,000 observed time units. No trajectories are combined.',fontsize=12,color=muted)
    fig.text(.04,.035,'The distinction is a row label and circulation sign. Symmetry-related outcomes are not separate discoveries.\nFinite sampled histories support this organization; they do not certify eight basins or an exact symmetry group.',fontsize=11,color=muted)
    fig.subplots_adjust(left=.01,right=.99,top=.855,bottom=.10,wspace=0,hspace=.08)
    fig.savefig(out/'eight-orientations.png',dpi=150,facecolor=bg);plt.close(fig)
    fig=plt.figure(figsize=(13,7),facecolor=bg)
    for k,(clouds,title) in enumerate([(positions,'Original row labels'),(aligned,'After the measured row relabeling')]):
        ax=fig.add_subplot(1,2,k+1,projection='3d',facecolor=bg)
        for cloud,color in zip(clouds,colors):ax.plot(*cloud.T,color=color,lw=.22,alpha=.40)
        style(ax,allpoints,27,-45,zoom=1.0);ax.set_title(title,fontsize=13)
    fig.suptitle('EIGHT SEPARATE HISTORIES | COMPARISON OVERLAY',x=.04,ha='left',y=.97,fontsize=19,weight='bold')
    fig.text(.04,.89,'The right panel only permutes the four physical row coordinates. It does not fit or deform the curves.',fontsize=11,color=muted)
    fig.text(.04,.045,'This overlay compares eight trajectories; neither panel represents a single trajectory.\nThe underlying comparison uses all c, s and X coordinates, not just this three-dimensional view.',fontsize=10,color=muted)
    fig.subplots_adjust(left=0,right=1,top=.82,bottom=.12,wspace=0)
    fig.savefig(out/'relabeling-overlay.png',dpi=150,facecolor=bg);plt.close(fig)
    write(out/'analysis.json',dict(declared=32,completed=32,classes=[dict(distinguished=k[0],circulation_sign=k[1],count=v) for k,v in sorted(counts.items())],
        best_median_range=[min(r['best']['median'] for r in runs),max(r['best']['median'] for r in runs)],
        best_p95_range=[min(r['best']['p95'] for r in runs),max(r['best']['p95'] for r in runs)],
        near_permutation_counts=dict(Counter(len(r['near_permutations']) for r in runs)),runs=runs,
        criterion='Near: median and p95 both at most three times that fresh start best-by-median permutation. Descriptive heuristic, not a proof.',
        comparison_direction='Queries from each fresh-start history to the fixed full reference cloud; not symmetric Hausdorff distance.',
        representatives='Lowest start ID in each (smallest mean position row, signed circulation of sorted remaining rows) class.',
        limitations='Finite histories and one initialization family. No basin-probability, attraction, exact isotropy, conjugacy or novelty proof.'))
    write(out/'render.json',dict(sources=hashes,contrast_matrix=H.tolist(),rule='Last600tau,stride2,uniform.1tau; shared bounds, camera27/-45, equal spatial axis units.',
        source_comparisons_sha256=hashlib.sha256((args.source/'comparisons.json').read_bytes()).hexdigest(),
        display_inputs_sha256=hashlib.sha256((out/'display-inputs.npz').read_bytes()).hexdigest()))
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print('classes',dict(counts),'near permutations',dict(Counter(len(r['near_permutations']) for r in runs)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);main(p.parse_args())
