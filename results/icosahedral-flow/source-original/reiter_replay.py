"""Independent implementation of Reiter's 1999 dodecahedral Figure 1 map.

Coefficients and checkpoints: author-hosted dodec02_pp.pdf, Appendices 1-2.
No source J implementation is copied; the mathematical formula is implemented.
"""
from pathlib import Path
import json,time
import numpy as np
from scipy.spatial.transform import Rotation
from numba import njit

root=Path(__file__).parent;out=root/'reiter-replay';out.mkdir(exist_ok=False)
coeff=np.array([
 .0180157305,.0903833937,-.1262988060,
 .0959720238,-.0981938805,-.0445451871,
 .1290451044,-.1389270186,-.0463975410,
 -.1034676072,.1199258382,-.1064360991,
 .0284841300,-.0672225822,-.0099387726,
 .0590534490,.1113198750,-.1468569510,
 -.1247709900,-.0260266215,-.0294227928,
 -.1088782995,.0824206119,.1432259172,
 -.0020062041,-.1182703734,.0298362312,
 -.1424609703,-.0415270104,-.1444622841,
 -.0776052450,.0886481373,.1092472107,
 .1178740806,-.0903261324,-.1113043062,
 .1085292588,.0512562201,-.1367073846,
 .0589886136,-.0783689961,-.1477134882,
 -.1205961711,-.0598447191,.0898105986,
 .1467316815,.1193712204,-.1278978093,
 -.0784800531,.0857490972,-.0149189145,
 .0578068920,-.1395622734,.0768726459,
 -.1014373995,.0416295027,.0670539840,
 -.1236866562,-.1016276886,.1434392721,
 -.0161524251,.0261961311,-.1216222695,
 .0945181746,.0669648237,-.1222031436,
 -.0682310184,.1412775246,-.0486418341,
 -.0233046243,.1191796536,-.0475573362,
 -.0961468965,-.1408879608,-.0039553878,
 .1218001257,-.1052849598,-.1243178283,
 .0902625387,-.0575116383,.0018951813]).reshape(3,3,3,3)
# Build the exact orientation used in the paper from its three generators.
phi=(1+np.sqrt(5))/2;inv=1/phi
generators=[.5*np.array([[1,phi,inv],[-phi,inv,1],[inv,-1,phi]]),.5*np.array([[phi,inv,-1],[inv,1,phi],[1,-phi,inv]]),np.diag([-1,1,1])]
group=[np.eye(3)];cursor=0
while cursor<len(group):
    g=group[cursor];cursor+=1
    for a in generators:
        b=g@a
        if min(np.linalg.norm(b-c) for c in group)>1e-10:group.append(b)
assert len(group)==120
group=np.array(group)
@njit
def step(x):
    result=np.zeros(3)
    for g in group:
        y=g@x;p=np.zeros(3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    p+=coeff[i,j,k]*y[0]**k*y[1]**j*y[2]**i
        result+=g.T@p
    return result

x=np.array([.1,.2,.3]);check=[x.copy()]
for _ in range(4):x=step(x);check.append(x.copy())
expected=np.array([[.1,.2,.3],[.253974,.507183,.761566],[.30382,.527163,.874296],[.234707,.414627,.717742],[.313538,.570422,.994186]])
error=float(np.max(np.abs(np.array(check)-expected)))
print(json.dumps({'checkpoint_error':error,'actual':np.array(check).tolist()}),flush=True)
assert error<6e-7,'Published checkpoints do not match; do not use the replay.'

@njit
def orbit(count):
    x=np.array([.1,.2,.3]);tail=np.empty((count,3))
    for i in range(count+10000):
        x=step(x)
        if i>=10000:tail[i-10000]=x
    return tail
begin=time.time();points=orbit(320000);assert np.isfinite(points).all()
np.savez_compressed(out/'trajectory.npz',points=points,group=group,coefficients=coeff)
(out/'checks.json').write_text(json.dumps({'source':'https://webbox.lafayette.edu/~reiterc/dodec/dodec02_pp.pdf','figure':1,'group_size':120,'checkpoint_error':error,'actual_checkpoints':np.array(check).tolist(),'warm':10000,'observed':320000,'seconds':time.time()-begin,'limits':'Only Figure 1 is reproduced. Figure 2 parameters were not supplied in this appendix. This is a geometric reference, not a validated chaos certificate.'},indent=2))
print('Replay complete',flush=True)
