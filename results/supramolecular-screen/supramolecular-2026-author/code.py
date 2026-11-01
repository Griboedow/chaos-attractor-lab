#!/usr/bin/env python
# coding: utf-8

import numpy as np

print("Assigning parameters")

# Simulation parameters
tmax = 100
nstep = 1000000
dt = tmax / nstep
dx = 0.5
x = np.arange(0, 100, dx)
space = len(x)
tim = nstep

# Model parameters
ka = 4
kd2 = 5000
kd3 = 0.75
kp = 1
ctot = 1e4
C =4.6

# Initial condition
print("Setting up initial condition")

# Steady state values
a1ss=5000
m1ss=86.0367
m0ss=1.96559

# (Random perturbation)
a1 = a1ss + a1ss* 1e-3 * np.random.rand(space, space)
m1 = m1ss + m1ss* 1e-3 * np.random.rand(space, space)
m0 = m0ss + m0ss* 1e-3 * np.random.rand(space, space)
d  = (ctot - a1 - m1) + (ctot-a1ss-m1ss)* 1e-3 * np.random.rand(space, space)

# (Gaussian perturbation)
a1 = a1ss + np.zeros((space, space))
xc = (space*dx)/2 
yc = (space*dx)/2 
inv_sigma = 0.5
for i in range(space): 
    for j in range (space):
        x = i*dx
        y = j*dx
        m1[i,j] = m1ss * 1e-3 * np.exp(- ((x-xc)**2 + (y-yc)**2)/2*inv_sigma**2) + m1ss
        m0[i,j] = m0ss * 1e-3 * np.exp(- ((x-xc)**2 + (y-yc)**2)/2*inv_sigma**2) + m0ss
d = (ctot-a1-m1) + np.zeros((space, space))

d1 = 1 # Monomer diffusion coefficient

# Operators 
def laplacian(Z):
    return (Z[0:-2,1:-1] + Z[2:,1:-1] + Z[1:-1,0:-2] + Z[1:-1,2:] - 4*Z[1:-1,1:-1]) / dx**2

def apply_bc(Z):
    Z[0, 1:-1] = Z[1, 1:-1]
    Z[-1, 1:-1] = Z[-2, 1:-1]
    Z[1:-1, 0] = Z[1:-1, 1]
    Z[1:-1, -1] = Z[1:-1, -2]
    Z[0,0], Z[0,-1], Z[-1,0], Z[-1,-1] = Z[1,1], Z[1,-2], Z[-2,1], Z[-2,-2]
    return Z

# Calculating RHS
def compute_rhs(d, a1, m1, m0):

    D = m0 / (m1 + 1e-12) # Polymer diffusion coefficient

    # Monomers 
    dd_rhs = (-ka*d[1:-1,1:-1] + 2*kd2*m0[1:-1,1:-1]) + d1*laplacian(d)

    aa1_rhs = (ka*d[1:-1,1:-1] - 2*kp*a1[1:-1,1:-1]*m0[1:-1,1:-1]) + d1*laplacian(a1)

    # Polymers
    lap_m1 = laplacian(m1)
    lap_m0 = laplacian(m0)

    grad_term_m1 = 0.25/dx**2 * ((D[2:,1:-1] - D[0:-2,1:-1]) * (m1[2:,1:-1] - m1[0:-2,1:-1]) +
        (D[1:-1,2:] - D[1:-1,0:-2]) * (m1[1:-1,2:] - m1[1:-1,0:-2]))

    grad_term_m0 = 0.25/dx**2 * ((D[2:,1:-1] - D[0:-2,1:-1]) * (m0[2:,1:-1] - m0[0:-2,1:-1]) +
        (D[1:-1,2:] - D[1:-1,0:-2]) * (m0[1:-1,2:] - m0[1:-1,0:-2]))

    mm1_rhs = (-2*kd2*m0[1:-1,1:-1] +
               2*kp*a1[1:-1,1:-1]*m0[1:-1,1:-1]) + \
              (D[1:-1,1:-1]*lap_m1 + grad_term_m1)

    mm0_rhs = (kd3*(m1[1:-1,1:-1] - 4*m0[1:-1,1:-1]) -
              ((2*kd2*C*m0[1:-1,1:-1]) /
               ((m1[1:-1,1:-1]/m0[1:-1,1:-1]) - 2)**2 +
               2*kd3*m0[1:-1,1:-1] +
               kp*m0[1:-1,1:-1]**2)) + \
              (D[1:-1,1:-1]*lap_m0 + grad_term_m0)

    return dd_rhs, aa1_rhs, mm1_rhs, mm0_rhs

# RK4 integration
print("Starting RK4 integration")
for k in range(1, tim):

    k1 = compute_rhs(d, a1, m1, m0)

    d2 = d.copy(); a12 = a1.copy(); m12 = m1.copy(); m02 = m0.copy()
    d2[1:-1,1:-1] += 0.5*dt*k1[0]
    a12[1:-1,1:-1] += 0.5*dt*k1[1]
    m12[1:-1,1:-1] += 0.5*dt*k1[2]
    m02[1:-1,1:-1] += 0.5*dt*k1[3]

    d2  = apply_bc(d2)
    a12 = apply_bc(a12)
    m12 = apply_bc(m12)
    m02 = apply_bc(m02)

    k2 = compute_rhs(d2, a12, m12, m02)

    d3 = d.copy(); a13 = a1.copy(); m13 = m1.copy(); m03 = m0.copy()
    d3[1:-1,1:-1] += 0.5*dt*k2[0]
    a13[1:-1,1:-1] += 0.5*dt*k2[1]
    m13[1:-1,1:-1] += 0.5*dt*k2[2]
    m03[1:-1,1:-1] += 0.5*dt*k2[3]

    d3  = apply_bc(d3)
    a13 = apply_bc(a13)
    m13 = apply_bc(m13)
    m03 = apply_bc(m03)

    k3 = compute_rhs(d3, a13, m13, m03)

    d4 = d.copy(); a14 = a1.copy(); m14 = m1.copy(); m04 = m0.copy()
    d4[1:-1,1:-1] += dt*k3[0]
    a14[1:-1,1:-1] += dt*k3[1]
    m14[1:-1,1:-1] += dt*k3[2]
    m04[1:-1,1:-1] += dt*k3[3]

    d4  = apply_bc(d4)
    a14 = apply_bc(a14)
    m14 = apply_bc(m14)
    m04 = apply_bc(m04)

    k4 = compute_rhs(d4, a14, m14, m04)

    d[1:-1,1:-1] += (dt/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
    a1[1:-1,1:-1] += (dt/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    m1[1:-1,1:-1] += (dt/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
    m0[1:-1,1:-1] += (dt/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3])

    d = apply_bc(d)
    a1 = apply_bc(a1)
    m1 = apply_bc(m1)
    m0 = apply_bc(m0)

    if k % 1000 == 0:
        print("Completed step", k)
        np.savez(f"step_{k/10000}.npz", d=d, a1=a1, m1=m1, m0=m0)

    d_mean.append(np.mean(d))
    a1_mean.append(np.mean(a1))
    m1_mean.append(np.mean(m1))
    m0_mean.append(np.mean(m0))

print("Completed integration")


