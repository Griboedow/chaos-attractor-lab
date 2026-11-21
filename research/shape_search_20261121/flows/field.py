"""Spatially periodic cross-couplings with diagonal damping in three-state space."""
import math

import numpy as np
from numba import njit


@njit(cache=True)
def velocity_tangent(state, p):
    result = np.empty(6)
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        d, a, b, c, phi, psi, chi, m, n = p[i]
        u = state[j] + phi
        v = state[k] + psi
        w = m * state[j] + n * state[k] + chi
        result[i] = -d * state[i] + a * math.sin(u) + b * math.sin(v) + c * math.sin(w)
        du = a * math.cos(u) + m * c * math.cos(w)
        dv = b * math.cos(v) + n * c * math.cos(w)
        result[i + 3] = -d * state[i + 3] + du * state[j + 3] + dv * state[k + 3]
    return result


@njit(cache=True, nogil=True)
def integrate(p, initial, dt, burn_steps, block_steps, blocks, stride, derivative=velocity_tangent):
    """Classical RK4, tangent rescaled every five steps and at block boundaries."""
    state = np.empty(6)
    state[:3] = initial
    state[3:] = np.array([1., 2., 3.]) / math.sqrt(14.)
    rates = np.empty(blocks)
    trace = np.empty((blocks * block_steps // stride, 3))
    sample = 0
    for block in range(-1, blocks):
        steps = burn_steps if block == -1 else block_steps
        growth = 0.
        for step in range(1, steps + 1):
            k1 = derivative(state, p)
            k2 = derivative(state + .5 * dt * k1, p)
            k3 = derivative(state + .5 * dt * k2, p)
            k4 = derivative(state + dt * k3, p)
            state += dt / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
            if step % 5 == 0 or step == steps:
                magnitude = math.sqrt(state[3]**2 + state[4]**2 + state[5]**2)
                growth += math.log(magnitude)
                state[3:] /= magnitude
            if block >= 0 and (block * block_steps + step) % stride == 0:
                trace[sample] = state[:3]
                sample += 1
        if block >= 0:
            rates[block] = growth / (steps * dt)
    return rates, trace


def make_case(seed, case):
    rng = np.random.default_rng(np.random.SeedSequence([seed, case]))
    variant = case % 6
    shared = [10**rng.uniform(-1.45, -.48), 1., rng.uniform(-1.25, 1.25),
              rng.uniform(-1.15, 1.15), 0., 0., 0., rng.choice([1, 2, 3]), rng.choice([-1, 1])]
    p = np.tile(shared, (3, 1))
    if variant in (1, 4):
        p[:, 4:7] = rng.uniform(-np.pi, np.pi, size=3)
    if variant in (2, 3, 5):
        p[:, 4:7] = rng.uniform(-np.pi, np.pi, size=(3, 3))
    if variant in (3, 4, 5):
        p[:, :4] *= rng.uniform(.65, 1.35, size=(3, 4))
        p[:, 7] = rng.choice([1, 2, 3], size=3)
        p[:, 8] = rng.choice([-1, 1], size=3)
    return p, rng.uniform(-4., 4., size=3)


@njit(cache=True)
def ring_tangent(state, p):
    """A saturated angular resonance with radial and axial feedback."""
    order = int(p[0])
    alpha, beta, omega, coupling, twist, strength, decay, radial, axial, phase = p[1:]
    u = complex(state[0], state[1])
    v = complex(state[3], state[4])
    z, vz = state[2], state[5]
    r2 = u.real*u.real + u.imag*u.imag
    dr2 = 2 * (u.real*v.real + u.imag*v.imag)
    factor = 2 / (1 + r2)
    factor_slope = -dr2 / (1 + r2)
    base = complex(alpha-r2+beta*z, omega+coupling*z+twist*r2)
    base_tangent = complex(-dr2+beta*vz, coupling*vz+twist*dr2)
    angular_scale = strength * factor**((order-2)/2)
    angular = u.conjugate()**(order-1)
    angular_tangent = (order-1) * u.conjugate()**(order-2) * v.conjugate()
    du = u * base + angular_scale * angular
    dv = v*base + u*base_tangent + angular_scale*(angular_tangent + angular*(order-2)/2*factor_slope)
    direction = complex(math.cos(phase), math.sin(phase))
    axial_scale = axial * factor**(order/2)
    drive = (direction*u**order).real
    drive_tangent = (direction * order * u**(order-1)*v).real
    dz = -decay*z + radial*(r2-1)/(1+r2) + axial_scale*drive
    dw = -decay*vz + 2*radial*dr2/(1+r2)**2 + axial_scale*(drive_tangent+drive*order/2*factor_slope)
    return np.array([du.real, du.imag, dz, dv.real, dv.imag, dw])


def make_ring(seed, case):
    rng = np.random.default_rng(np.random.SeedSequence([seed, case]))
    order = 3 + case % 6
    p = np.array([order, rng.uniform(.4, 1.4), rng.uniform(-1.6, 1.6), rng.uniform(-2, 2),
                  rng.uniform(-3., 3.), rng.uniform(-2., 2.), rng.uniform(.12, 1.2),
                  10**rng.uniform(-.8, .25), rng.uniform(-2., 2.), rng.uniform(-2., 2.),
                  rng.uniform(-np.pi, np.pi)])
    angle = rng.uniform(0, 2*np.pi)
    initial = np.array([np.cos(angle), np.sin(angle), rng.uniform(-.8, .8)])
    return p, initial


@njit(cache=True)
def rotor_tangent(state, p):
    """Damped torque-driven precession with a quartic angular potential."""
    x, y, z = state[:3]
    u, v, w = state[3:]
    d = p[:3]
    torque = p[3:6]
    a = p[6:9]
    q = p[9:12]
    h = p[12:15]
    c = p[15]
    gx = a[0]+q[0]*x+h[0]*x**3+c*y*z
    gy = a[1]+q[1]*y+h[1]*y**3+c*x*z
    gz = a[2]+q[2]*z+h[2]*z**3+c*x*y
    dgx = (q[0]+3*h[0]*x*x)*u+c*(v*z+y*w)
    dgy = (q[1]+3*h[1]*y*y)*v+c*(u*z+x*w)
    dgz = (q[2]+3*h[2]*z*z)*w+c*(u*y+x*v)
    fx, fy, fz = y*gz-z*gy, z*gx-x*gz, x*gy-y*gx
    factor = 1/(1+x*x+y*y+z*z)
    slope = -2*(x*u+y*v+z*w)*factor**2
    return np.array([torque[0]-d[0]*x+factor*fx,
                     torque[1]-d[1]*y+factor*fy,
                     torque[2]-d[2]*z+factor*fz,
                     -d[0]*u+factor*(v*gz+y*dgz-w*gy-z*dgy)+slope*fx,
                     -d[1]*v+factor*(w*gx+z*dgx-u*gz-x*dgz)+slope*fy,
                     -d[2]*w+factor*(u*gy+x*dgy-v*gx-y*dgx)+slope*fz])


def make_rotor(seed, case):
    rng = np.random.default_rng(np.random.SeedSequence([seed, case]))
    d = np.full(3, 10**rng.uniform(-1.65,-.35)) * rng.uniform(.55,1.45,3)
    torque = rng.normal(size=3)
    torque *= d.mean()*rng.uniform(1.2,5)/np.linalg.norm(torque)
    p = np.r_[d, torque, rng.uniform(-2.,2.,3), rng.uniform(-4.,4.,3),
              rng.uniform(-3.,3.,3), rng.uniform(-8.,8.)]
    return p, rng.uniform(-2., 2., 3)


@njit(cache=True)
def plane_tangent(state, p):
    """A leaky three-state feedback network of smooth switching planes."""
    count = int(p[0])
    out = np.empty(6)
    for i in range(3):
        out[i] = p[4+i]-p[1+i]*state[i]
        out[3+i] = -p[1+i]*state[i+3]
    for r in range(count):
        j = 7 + 7*r
        value = math.tanh(p[j]*state[0]+p[j+1]*state[1]+p[j+2]*state[2]+p[j+3])
        diff = (1-value*value)*(p[j]*state[3]+p[j+1]*state[4]+p[j+2]*state[5])
        for i in range(3):
            out[i] += p[j+4+i]*value
            out[i+3] += p[j+4+i]*diff
    return out


def make_planes(seed, case):
    rng = np.random.default_rng(np.random.SeedSequence([seed, case]))
    count = int(rng.choice([5,7,9,11]))
    d = 10**rng.uniform(-1.3,-.15,size=3)
    normals = rng.normal(size=(count,3))
    normals *= (10**rng.uniform(-.35,.65,count)/np.linalg.norm(normals,axis=1))[:,None]
    bias = np.zeros(count) if case%3==0 else rng.uniform(-1.7,1.7,count)
    weights = rng.normal(size=(3,count)) * 10**rng.uniform(-.3,.6)/np.sqrt(count)
    basis, _ = np.linalg.qr(rng.normal(size=(3,3)))
    axial = 10**rng.uniform(-.8,.35)
    transverse = 10**rng.uniform(-1.,.12)
    frequency = 10**rng.uniform(-.12,.55)
    if case%2:
        axial *= -1
    else:
        transverse *= -1
    block = np.array([[axial,0.,0.],[0.,transverse,-frequency],[0.,frequency,transverse]])
    target = np.diag(d)+basis@block@basis.T
    slopes = (1-np.tanh(bias)**2)[:,None]*normals
    weights += (target-weights@slopes)@np.linalg.pinv(slopes)
    drift = -(weights @ np.tanh(bias))
    rows = np.column_stack((normals,bias,weights.T))
    return np.r_[count,d,drift,rows.ravel()], rng.uniform(-2.,2.,size=3)


@njit(cache=True)
def crossed_tangent(state, p):
    """Add resonant vector fields about two skew axes in the same XYZ space."""
    rotation = p[22:31].reshape(3,3)
    center = p[31:34]
    weight = p[34]
    local = np.empty(6)
    local[:3] = rotation.T @ (state[:3]-center)
    local[3:] = rotation.T @ state[3:]
    a = ring_tangent(state,p[:11])
    b = ring_tangent(local,p[11:22])
    a[:3] = (1-weight)*a[:3]+weight*(rotation @ b[:3])
    a[3:] = (1-weight)*a[3:]+weight*(rotation @ b[3:])
    return a


def make_crossed(seed, case):
    rng = np.random.default_rng(np.random.SeedSequence([seed,case]))
    base = [21253,20897,20738,21730,21001,20416,24830,24598,27124,25718,23772,23870]
    id_a,id_b = rng.choice(base,2,replace=False)
    p_a,start = make_ring(seed,int(id_a))
    p_b,_ = make_ring(seed,int(id_b))
    theta,psi = rng.uniform(0,2*np.pi,2)
    tilt = rng.uniform(.12,1.5)
    def zrot(angle):
        return np.array([[np.cos(angle),-np.sin(angle),0],
                         [np.sin(angle),np.cos(angle),0],[0,0,1]])
    yrot = np.array([[np.cos(tilt),0,np.sin(tilt)],[0,1,0],[-np.sin(tilt),0,np.cos(tilt)]])
    rotation = zrot(theta) @ yrot @ zrot(psi)
    center = rng.uniform(-.5,.5,size=3)
    weight = 10**rng.uniform(-1.4,-.21)
    return np.r_[p_a,p_b,rotation.ravel(),center,weight],start
