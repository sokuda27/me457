# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 10:57:19 2025

@author: saira
"""

import numpy as np
import matplotlib.pyplot as plt
import integrators as intg

# State Initial Conditions                     
p_n0 = 1; p_e0 = 0; p_d0 = 2
u_0 = 1; v_0 = 1; w_0 = 0
phi_0 = 0; theta_0 = 0; psi_0 = 0
p_0 = 0; q_0 = 0; r_0 = 0

# System conditions
# wind in body frame
u_w0 = 10
v_w0 = 10
w_w0 = 0

# Input Initial Conditions (forces and moments)
#Fx_0 = 1; Fy_0 = 1; Fz_0 = 1
#lx_0 = 1; my_0 = 1; nz_0 = 1

#physical parameters
mass = 11
g = 9.8
rho = 1.268
Jx=.824; Jy=1.135; Jz=1.759
Jxy=0; Jyz=0; Jxz=0
S = .55; b = 2.9; c = .19; e = .9
omega = 44.4/.0659
D = .508

T = Jx*Jz - Jxz**2
T1 = Jxz*(Jx-Jy+Jz)/T
T2 = (Jz*(Jz-Jy) + Jxz**2)/T
T3 = Jz/T
T4 = Jxz/T
T5 = (Jz-Jx)/Jy
T6 = Jxz/Jy
T7 = ((Jx-Jy)*Jx + Jxz**2)/T
T8 = Jx/T

# Longitudinal coeffs
C_L_0 = .23; C_D_0 = .043; C_m_0 = .0135; C_L_a = 5.61; C_D_a = .03; C_m_a = -2.74; C_L_q = 7.95; C_D_q = 0
C_m_q = -38.21; C_L_pe = .13; C_D_pe = .0135; C_m_pe = -.99; 
# Laterial Coeffs
C_Y_0 = 0; C_l_0 = 0; C_n_0 = 0; C_Y_B = -.83; C_l_B = -.13; C_n_B = .073; C_Y_p = 0; C_l_p = -.51; C_n_p = -.069
C_Y_r = 0; C_l_r = .25; C_n_r = -.095; C_Y_pa = .075; C_l_pa = .17; C_n_pa = -.011; C_Y_pr = .19; C_l_pr = .0024; C_n_pr = -.069
#Aero coeffs
CT2 = -0.1079; CT1 = -0.06044; CT0 = 0.09357
CQ1 = 0.004970; CQ2 = -0.01664; CQ0 = 0.005230

# Nonlinear state space form:
#  xdot = f(t, x, u), 
#   t: time
#   x: state vector
#   u: input

def f(t, x, u_i):
    p_n, p_e, p_d, u, v, w, phi, theta, psi, p, q, r = x #state variables
    u_w, v_w, w_w = u_i #input variables
    
    
    # variables
    u_r = u-u_w; v_r = v-v_w; w_r = w-w_w
    #print(u_w)
    
    V_a = np.sqrt(u_r**2 + v_r**2 + w_r**2)
    alpha = np.arctan(w_r/u_r)
    beta = np.arcsin(v_r/V_a)
    pa = p; pe = -q; pr = -r #partial a,e,r  corresponding to aileron, elevator, and rudder
    J = V_a*2*np.pi/(omega*D)
    
    #linear lift and drag models
    C_L = C_L_0 + C_L_a*alpha
    C_D = C_D_0 + C_D_a*alpha
    
    #coefficients of drag
    C_x_alpha = -C_D*np.cos(alpha) + C_L*np.sin(alpha)
    C_x_q_alpha = -C_D_q*np.cos(alpha) + C_L_q*np.sin(alpha)
    C_x_p_e_alpha = -C_D_pe*np.cos(alpha) + C_L_pe*np.sin(alpha)
    C_z_alpha = -C_D*np.sin(alpha) - C_L*np.cos(alpha)
    C_z_q_alpha = -C_D_q*np.sin(alpha) - C_L_q*np.cos(alpha)
    C_z_p_e_alpha = -C_D_pe*np.sin(alpha) - C_L_pe*np.cos(alpha)
    
    #aerodynamic coefficients
    C_T_J = CT2*J**2 + CT1*J + CT0
    C_Q_J = CQ2*J**2 + CQ1*J + CQ0
    
    #propeller thrust and torque
    T_p = rho*(omega/(2*np.pi))**2*(D**4)*C_T_J
    Q_p = rho*(omega/(2*np.pi))**2*(D**5)*C_Q_J
    
    #new forces and moments with wind and drag and propulsion
    fx = -mass*g*np.sin(theta) + (.5*rho*V_a**2*S)*(C_x_alpha + C_x_q_alpha*(c/(2*V_a))*q + C_x_p_e_alpha*pe) + T_p
    fy = mass*g*np.cos(theta)*np.sin(phi) + (.5*rho*V_a**2*S)*(C_Y_0 + C_Y_B + C_Y_p*(b/(2*V_a))*p + C_Y_r*(b/(2*V_a))*r + C_Y_pa*pa + C_Y_pr*pr) + 0
    fz = mass*g*np.cos(theta)*np.cos(phi) + (.5*rho*V_a**2*S)*(C_z_alpha + C_z_q_alpha*(c/(2*V_a))*q + C_z_p_e_alpha*pe) + 0
    l = (.5*rho*V_a**2*S)*b*(C_l_0 + C_l_B*beta + C_l_p*(b/(2*V_a))*p + C_l_r*(b/(2*V_a))*r + C_l_pa*pa + C_l_pr*pr) + Q_p
    m = (.5*rho*V_a**2*S)*c*(C_m_0 + C_m_a*alpha + C_m_q*(c/(2*V_a))*q + C_m_pe*pe) + 0
    n = (.5*rho*V_a**2*S)*b*(C_n_0 + C_n_B*beta + C_n_p*(b/(2*V_a))*p + C_n_r*(b/(2*V_a))*r + C_n_pa*pa + C_n_pr*pr) + 0
    
    # equations of motion
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
    #print(psi)
    return np.array([p_ndot, p_edot, p_ddot,udot, vdot, wdot, phidot, thetadot, psidot, pdot, qdot, rdot])


t = 0
x = np.array([p_n0, p_e0, p_d0, u_0, v_0, w_0, phi_0, theta_0, psi_0, p_0, q_0, r_0])
#u_i = np.array([Fx_0, Fy_0, Fz_0, lx_0, my_0, nz_0])
u_i = np.array([u_w0, v_w0, w_w0])
dt = .01; n = 10
print(u_i)

integrator = intg.RungeKutta4(dt, f)

t_history = [0]
x_history = [x]
for i in range(n):
    
    x = integrator.step(t, x, u_i)
    t = (i+1) * dt
    t_history.append(t)
    x_history.append(x)

#print(x_history)
intg.__doc__
#plt.figure()
#plt.title('pn')
p_north = [arr[0] for arr in x_history]
p_east = [arr[1] for arr in x_history]
p_down = [arr[2] for arr in x_history]
#plt.plot(t_history, p_north)
#plt.legend(["pn", "pe", "pd"])
#plt.show()

plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(p_north, p_east, p_down)
ax.set_title('Path of aircraft')
ax.set_xlabel('north', fontsize=10)
ax.set_ylabel('east', fontsize=10)
ax.set_zlabel('down', fontsize=10)
plt.show()

