"""Numerical multiple-shooting realizations of binary covering itineraries."""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import root
from continue_raar_cycles import iterate

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);data=json.loads(args.input.read_text());centers=np.array(data['centers']);S=np.array(data['bases'])*np.array(data['radii'])[:,None,:];C=np.linalg.inv(S);p=np.array(data['params']);rows=[];begin=time.time()
    for word in ['1','01','001','011','0001','0011','001011']:
        labels=np.array([v for letter in word for v in ([0]*7 if letter=='0' else list(range(7)))]);n=len(labels)
        def equations(flat):
            q=flat.reshape(n,3);f=np.empty_like(q);jac=np.zeros((3*n,3*n))
            for i,label in enumerate(labels):
                nxt=(i+1)%n;target=labels[nxt];x=centers[label]+S[label]@q[i];y,J,_,_,status=iterate(x,p,7);assert status==0;f[i]=C[target]@(y-centers[target])-q[nxt];jac[3*i:3*i+3,3*i:3*i+3]=C[target]@J@S[label];jac[3*i:3*i+3,3*nxt:3*nxt+3]-=np.eye(3)
            return f.ravel(),jac
        solution=root(equations,np.zeros(3*n),jac=True,method='hybr',options={'xtol':1e-10,'maxfev':400});q=solution.x.reshape(n,3);states=np.array([centers[label]+S[label]@v for label,v in zip(labels,q)]);orbit=[];monodromy=np.eye(3);residuals=[];axes=[]
        for i,x in enumerate(states):
            y,M,gap,segment,status=iterate(x,p,7);orbit.extend(segment);monodromy=M@monodromy;residuals.append(np.max(abs(y-states[(i+1)%n])));axes.append(gap)
        maximum_coordinate=float(np.max(abs(q)));maximum_residual=float(max(residuals));good=maximum_coordinate<1 and maximum_residual<1e-10;multipliers=np.linalg.eigvals(monodromy)
        np.savez_compressed(args.out/f'word-{word}.npz',params=p,word=word,labels=labels,coordinates=q,return_states=states,orbit=np.array(orbit),residuals=residuals,multipliers=multipliers)
        row=dict(word=word,steps=49*len(word),solver_success=bool(solution.success),solver_message=solution.message,maximum_window_coordinate=maximum_coordinate,maximum_physical_residual=maximum_residual,minimum_axis=float(min(axes)),accepted_numerical=good,multipliers=[[float(v.real),float(v.imag)] for v in multipliers]);rows.append(row);print(json.dumps(row),flush=True)
        if word=='1' and good:(args.out/'certificate-input.json').write_text(json.dumps(dict(cases=[dict(label='symbolic-period-49',period=49,params=p.tolist(),state=states[0].tolist())]),indent=2)+'\n',encoding='utf-8')
    (args.out/'report.json').write_text(json.dumps(dict(cases=rows,seconds=time.time()-begin,scope='Numerical multiple-shooting orbits inside the specified windows. Periods and stability require separate certification; failed solver flags are retained even when physical residuals are small.'),indent=2)+'\n',encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
