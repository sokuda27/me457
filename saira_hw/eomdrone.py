# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 10:57:19 2025

@author: saira
"""

import numpy as np
import matplotlib.pyplot as plt
import integrators as intg


# Nonlinear state space form:
#  xdot = f(t, x, u), 
#   t: time
#   x: state vector
#   u: input

def f(t, x, u_i):
    p_n, p_e, p_d, u, v, w, phi, theta, psi, p, q, r = x #state variables
    fx, fy, fz, l, m, n = u_i #input variables
    
    p_ndot = (np.cos(theta)*np.cos(psi))*u + (np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi))*v + (np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi))*w
    p_edot = (np.cos(theta)*np.sin(psi))*u + (np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi))*v + (np.cos(phi)*np.sin(theta)*np.sin(psi) - np.sin(phi)*np.cos(psi))*w     
    p_ddot = (-np.sin(theta))*u + (np.sin(phi)*np.cos(theta))*v + (np.cos(phi)*np.cos(theta))*w
    udot = r*v - q*w + (1/mass)*fx 
    vdot = p*w - r*u + (1/mass)*fy 
    wdot = q*u - p*v + (1/mass)*fz 
    phidot = 1*p + np.sin(phi)*np.tan(theta)*q + np.cos(phi)*np.tan(theta)*r 
    thetadot = np.cos(phi)*q - np.sin(phi)*r 
    psidot = (np.sin(phi)/np.cos(theta))*q + (np.cos(phi)/np.cos(theta))*r 
    pdot = (T1*p*q -T2*q*r) + (T3*l + T4*n)
    qdot = (T5*p*r - T6*(p**2 - r**2)) + (1/Jy)*m
    rdot = (T7*p*q - T1*q*r) + (T4*l + T8*n)
    
    return np.array([p_ndot, p_edot, p_ddot,udot, vdot, wdot, phidot, thetadot, psidot, pdot, qdot, rdot])

# State Initial Conditions                     
p_n0 = 0; p_e0 = 0; p_d0 = 100
u_0 = 10; v_0 = 0; w_0 = 0
phi_0 = 0; theta_0 = 0; psi_0 = 0
p_0 = 0; q_0 = 0; r_0 = 0

# Input Initial Conditions (forces and moments)
Fx_0 = 1; Fy_0 = 0; Fz_0 = 9.8;
lx_0 = 0; my_0 = 0; nz_0 = 0

#moments of inertia
mass = 1
Jx=1; Jy=1; Jz=2
Jxy=0; Jyz=0; Jxz=0

T = Jx*Jz - Jxz**2
T1 = Jxz*(Jx-Jy+Jz)/T
T2 = (Jz*(Jz-Jy) + Jxz**2)/T
T3 = Jz/T
T4 = Jxz/T
T5 = (Jz-Jx)/Jy
T6 = Jxz/Jy
T7 = ((Jx-Jy)*Jx + Jxz**2)/T
T8 = Jx/T


t = 0
x = np.array([p_n0, p_e0, p_d0, u_0, v_0, w_0, phi_0, theta_0, psi_0, p_0, q_0, r_0])
u_i = np.array([Fx_0, Fy_0, Fz_0, lx_0, my_0, nz_0])
dt = .1; n = 100


integrator = intg.RungeKutta4(dt, f)

t_history = [0]
x_history = [x]
for i in range(n):
    
    x = integrator.step(t, x, u_i)
    t = (i+1) * dt
    t_history.append(t)
    x_history.append(x)


intg.__doc__
plt.figure()
plt.title('pn')
p_north = [arr[0] for arr in x_history]
p_east = [arr[1] for arr in x_history]
p_down = [arr[2] for arr in x_history]
plt.plot(t_history, p_north)
#plt.legend(["pn", "pe", "pd"])
plt.show()

plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(p_north, p_east, p_down)
ax.set_title('Path of aircraft')
ax.set_xlabel('north', fontsize=10)
ax.set_ylabel('east', fontsize=10)
ax.set_zlabel('down', fontsize=10)
plt.show()
