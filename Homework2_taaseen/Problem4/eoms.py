# In[ ]:
# Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import integrators as intg

######################################################################
# Author: Taaseen Jahan                                              #
#                                                                    #
# Implements the equations of motion for a rigid body drone in 3D.   #
# The rigid body is assumed to be the Aerosonde with parameters      #
# taken from Small Unmanned Aircraft: Theory and Practice Supplement #
# by Beard and McLain.                                               #
#                                                                    #
# Received help and adapted code from Saira Billah.                  #
######################################################################

# Inputs
## Forces
f_x = 10 # Force in x
f_y = 10 # Force in y
f_z = 10 # Force in z

## Moments
l = 0.1  # Moment about x
m = 0.1  # Moment about y
n = 0.1  # Moment about z

## Input Arrays
f_b = np.array([f_x, f_y, f_z]) # Force Array
m_b = np.array([l, m, n])       # Moment Array

# Parameters
M = 1          # Mass in kg
J_x = 1        # Moment of Inertia about x in kg*m^2
J_y = 1        # Moment of Inertia about y in kg*m^2
J_z = 2        # Moment of Inertia about z in kg*m^2
J_xy = 0       # Product of Inertia in kg*m^2 (Assumed to be 0 for symmetric body)
J_yz = 0       # Product of Inertia in kg*m^2 (Assumed to be 0 for symmetric body)
J_xz = 0       # Product of Inertia in kg*m^2

J = np.array([[J_x, -J_xy, -J_xz],  # Inertia Matrix
              [-J_xy, J_y, -J_yz],
              [-J_xz, -J_yz, J_z]])

# Initial Conditions
## Translational Kinematics State
p_n_0 = 1
p_e_0 = 1
p_d_0 = 1
pos_0 = np.array([p_n_0, p_e_0, p_d_0]) # Position Vector
## Rotational Kinematics State
psi_0 = np.deg2rad(1)
theta_0 = np.deg2rad(1)
phi_0 = np.deg2rad(1)
euler_0 = np.array([psi_0, theta_0, phi_0]) # Euler Angles
## Translational Dynamics State
u_0 = 1
v_0 = 1
w_0 = 1
V_0 = np.array([u_0, v_0, w_0]) # Linear Velocity Vector
## Rotational Dynamics State
p_0 = 1
q_0 = 0
r_0 = 2
W_0 = np.array([p_0, q_0, r_0]) # Angular Velocity Vector

IC = np.array([p_n_0, p_e_0, p_d_0, psi_0, theta_0, phi_0, u_0, v_0, w_0, p_0, q_0, r_0]) # Initial Conditions Vector

def f(t, x, u_i):
    p_n, p_e, p_d, psi, theta, phi, u, v, w, p, q, r = x
    
    # Rotation Matrices
    pos_rot_matrix = np.array([[np.cos(theta)*np.cos(psi), np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi)],
                             [np.cos(theta)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.cos(psi) + np.cos(phi)*np.cos(psi), np.cos(phi)*np.sin(theta)*np.sin(psi) + np.sin(phi)*np.cos(psi)],
                             [-np.sin(theta), np.sin(phi)*np.cos(theta), np.cos(phi)*np.cos(theta)]])
    euler_rot_matrix = np.array([[1, np.sin(phi)*np.tan(theta), np.cos(phi)*np.tan(theta)],
                                [0, np.cos(phi), -np.sin(phi)],
                                [0, np.sin(phi)/np.cos(theta), np.cos(phi)/np.cos(theta)]])
    
    # Gamma Values
    Gamma = J_x*J_z-J_xz**2
    Gamma1 = J_xz*(J_x-J_y+J_z)/Gamma
    Gamma2 = (J_z*(J_z-J_y)+J_xz**2)/Gamma
    Gamma3 = J_z/Gamma
    Gamma4 = J_xz/Gamma
    Gamma5 = (J_z-J_x)/J_y
    Gamma6 = J_xz/J_y
    Gamma7 = ((J_x-J_y)*J_x+J_xz**2)/Gamma
    Gamma8 = J_x/Gamma
    
    pos_dot = np.dot(pos_rot_matrix, np.array([u, v, w])) # From Slide 15
    V_dot = np.array([r*v-q*w, p*w-r*u, q*u-p*v]) + (1/M)*f_b # From Slide 15
    euler_dot = np.dot(euler_rot_matrix, np.array([p, q, r])) # From Slide 15
    W_dot = np.array([Gamma1*p*q-Gamma2*q*r, Gamma5*p*r-Gamma6*(p**2-r**2), Gamma7*p*q-Gamma1*q*r]) + np.array([Gamma3*l+Gamma4*n, (1/J_y)*M, Gamma4*l+Gamma8*n]) # From Slide 15
    
    return np.concatenate((pos_dot, V_dot, euler_dot, W_dot))
# Time Stepping Parameters
dt = 0.1; n = 10

t = 0; x = IC; u_i = 0
integrator = intg.RK4(dt, f)
t_history = [0]
x_history = [x]
for i in range(n):
    
    x = integrator.step(t, x, u_i)
    t = (i+1) * dt

    t_history.append(t)
    x_history.append(x)

intg.__doc__
plt.figure()
plt.title('Graph')
plt.plot(t_history, x_history)
plt.show()
