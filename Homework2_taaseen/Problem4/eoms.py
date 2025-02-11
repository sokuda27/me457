# Imports
import numpy as np
import matplotlib.pyplot as plt

######################################################################
# Author: Taaseen Jahan                                              #
#                                                                    #
# Implements the equations of motion for a rigid body drone in 3D.   #
# The rigid body is assumed to be the Aerosonde with parameters      #
# taken from Small Unmanned Aircraft: Theory and Practice Supplement #
# by Beard and McLain.                                               #
#                                                                    #
######################################################################

# Inputs
## Forces
f_x = 10 # Force in x
f_y = 10 # Force in y
f_z = 10 # Force in z

## Moments
l = 10   # Moment about x
m = 10   # Moment about y
n = 10   # Moment about z

## Input Arrays
f_b = np.array([f_x, f_y, f_z]) # Force Array
m_b = np.array([l, m, n])       # Moment Array

# Parameters
M = 11             # Mass in kg
J_x = 0.824        # Moment of Inertia about x in kg*m^2
J_y = 1.135        # Moment of Inertia about y in kg*m^2
J_z = 1.759        # Moment of Inertia about z in kg*m^2
J_xy = 0           # Product of Inertia in kg*m^2 (Assumed to be 0 for symmetric body)
J_yz = 0           # Product of Inertia in kg*m^2 (Assumed to be 0 for symmetric body)
J_xz = 0.120       # Product of Inertia in kg*m^2

J = np.array([[J_x, -J_xy, -J_xz],  # Inertia Matrix
              [-J_xy, J_y, -J_yz],
              [-J_xz, -J_yz, J_z]])

# Initial Conditions
## Translational Kinematics State
p_n = 0
p_e = 0
p_d = 0
## Rotational Kinematics State
psi = 0
theta = 0
phi = 0
## Translational Dynamics State
u = 0
v = 0
w = 0
## Rotational Dynamics State
p = 0
q = 0
r = 0
