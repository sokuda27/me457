# In[ ]:
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
u = 1
v = 1
w = 1
p = 1
q = 1
r = 1
V = np.array([u, v, w]) # Linear Velocity Vector
w = np.array([p, q, r]) # Angular Velocity Vector
w_bb = np.array([[0, -r, q],
                 [r, 0, -p],
                 [-q, p, 0]])

# Solve EOMs
V_solution = (1/M**2)*np.dot(-w_bb, V) + f_b
w_solution = np.dot(np.linalg.inv(J), (m_b-np.dot(w_bb,np.dot(J,w))))

# Print results
print('udot, vdot, wdot:')
print(V_solution)
print('pdot, qdot, rdot')
print(w_solution)
