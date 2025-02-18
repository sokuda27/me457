# In[ ]:
# Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

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
p_n = 1
p_e = 1
p_d = 1
pos = np.array([p_n, p_e, p_d]) # Position Vector
## Rotational Kinematics State
psi = np.deg2rad(1)
theta = np.deg2rad(1)
phi = np.deg2rad(1)
euler = np.array([psi, theta, phi]) # Euler Angles
## Translational Dynamics State
u = 1
v = 1
w = 1
V = np.array([u, v, w]) # Linear Velocity Vector
## Rotational Dynamics State
p = 1
q = 0
r = 2
w = np.array([p, q, r]) # Angular Velocity Vector

w_bb = np.array([[0, -r, q],
                 [r, 0, -p],
                 [-q, p, 0]])

# Solve EOMs
pos_rot_matrix = np.array([[np.cos(theta)*np.cos(psi), np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi)],
                         [np.cos(theta)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.cos(psi) + np.cos(phi)*np.cos(psi), np.cos(phi)*np.sin(theta)*np.sin(psi) + np.sin(phi)*np.cos(psi)],
                         [-np.sin(theta), np.sin(phi)*np.cos(theta), np.cos(phi)*np.cos(theta)]])
# pos_solution = np.dot(pos_rot_matrix, pos)
# V_solution = (1/M**2)*np.dot(-w_bb, V) + f_b
euler_rot_matrix = np.array([[1, np.sin(phi)*np.tan(theta), np.cos(phi)*np.tan(theta)],
                            [0, np.cos(phi), -np.sin(phi)],
                            [0, np.sin(phi)/np.cos(theta), np.cos(phi)/np.cos(theta)]])
# euler_solution = np.dot(euler_rot_matrix, euler)
# w_solution = np.dot(np.linalg.inv(J), (m_b-np.dot(w_bb,np.dot(J,w))))

pos_dot = pos_rot_matrix*pos
V_dot = (1/M)*(np.cross(-w,V) + f_b)
euler_dot = euler_rot_matrix*euler
w_dot= np.linalg.inv(J)*(M-np.cross(w, np.cross(J,w)))

dots = np.concatenate(pos_dot, V_dot, euler_dot, w_dot)
params = np.concatenate(pos, V, euler, w)

def model(y, t, params):
    """
    Defines the system of 12 ODEs.

    Args:
        y (array): A state vector with 12 elements.
        t (float): Time.
        params (tuple): A tuple of parameters for the model.

    Returns:
        dydt (array): An array of derivatives for each state variable.
    """
    # Example: Linear system with interaction
    A = params
    dydt = np.dot(A, y)
    return dydt

solve = odeint(model, y0, t, args=(dots,))