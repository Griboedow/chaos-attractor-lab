"""Generate an algebraically equivalent drag-Jacobian assembly for comparison."""
from pathlib import Path
root=Path(__file__).resolve().parent
source=(root/'spatial_fast.py').read_text(encoding='utf-8')
begin=source.index('@njit(cache=True)\ndef drag_jacobian(')
end=source.index('@njit(cache=True)\ndef hinge_eval(',begin)
replacement='''@njit(cache=True)
def drag_jacobian(q,nodes_old,nodes_prev,zeta,xi,dt,order):
    # The translational mass is symmetric. Precontract each tail with its
    # right tangent action, then fill only its two active rotational axes.
    # Drag-coefficient variations depend on the varied segment, so their
    # force vectors are formed once before propagation to upstream segments.
    frames,nodes=geometry(q)
    n=len(q);h=1.0/n;factor=(1.0 if order==1 else 1.5)/dt
    velocity=(nodes-nodes_old)/dt if order==1 else (1.5*nodes-2*nodes_old+0.5*nodes_prev)/dt
    ez=np.array([0.0,0.0,1.0]);s=-skew(ez)
    drags=np.empty((n,3,3));actions=np.empty((n,3,3));tails=np.zeros((n+1,3,3));forces=np.zeros((n+1,3))
    for i in range(n-1,-1,-1):
        drags[i]=(frames[i]*zeta)@frames[i].T*xi
        actions[i]=frames[i]@s
        tails[i]=tails[i+1]+drags[i]
        forces[i]=forces[i+1]+h*drags[i]@((velocity[i]+velocity[i+1])/2)
    jac=np.zeros((3*n,3*n))
    for j in range(n):
        right=(tails[j+1]+drags[j]/2)@actions[j]
        diagonal=(tails[j+1]+drags[j]/3)@actions[j]
        for i in range(j+1):
            block=diagonal if i==j else right
            for a in range(2):
                for b in range(2):
                    value=0.0
                    for component in range(3):
                        value+=actions[i,component,a]*block[component,b]
                    value*=factor*h**3
                    jac[3*i+a,3*j+b]=value
                    jac[3*j+b,3*i+a]=value
        weighted=velocity[j]/6+velocity[j+1]/3
        body=frames[j].T@(forces[j+1]+h*drags[j]@weighted)
        jac[3*j,3*j]-=h*body[2]
        jac[3*j+1,3*j+1]-=h*body[2]
        jac[3*j,3*j+2]+=h*body[0]
        jac[3*j+1,3*j+2]+=h*body[1]
        for axis in range(3):
            basis=np.zeros(3);basis[axis]=1.0
            cross=skew(basis)
            local=cross*np.expand_dims(zeta,0)-np.expand_dims(zeta,1)*cross
            delta=xi*frames[j]@local@frames[j].T
            avg_force=delta@((velocity[j]+velocity[j+1])/2)
            weighted_force=delta@weighted
            for i in range(j+1):
                force=weighted_force if i==j else avg_force
                for a in range(2):
                    value=0.0
                    for component in range(3):
                        value+=actions[i,component,a]*force[component]
                    jac[3*i+a,3*j+axis]+=h*h*value
    return jac


'''
(root/'spatial_accelerated.py').write_text(source[:begin]+replacement+source[end:],encoding='utf-8')
driver=(root/'restart.py').read_text(encoding='utf-8').replace('import spatial_fast as m','import spatial_accelerated as m')
driver=driver.replace("pars=report['parameters']", "pars=report['parameters'].copy()\n    if args.chi is not None:pars['chi']=args.chi")
driver=driver.replace("segments=len(q0),time=origin)","segments=len(q0),time=origin,parameters=report['parameters'])")
driver=driver.replace("ap.add_argument('--tolerance'", "ap.add_argument('--chi',type=float)\n    ap.add_argument('--tolerance'")
(root/'restart_accelerated.py').write_text(driver,encoding='utf-8')
