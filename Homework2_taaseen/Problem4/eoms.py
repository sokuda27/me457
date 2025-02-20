# In[ ]:
# Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import integrators as intg
from mpl_toolkits.mplot3d import Axes3D

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
f_x_0 = 1 # Force in x
f_y_0 = 0 # Force in y
f_z_0 = 9.81 # Force in z

## Moments
l_0 = 0  # Moment about x
m_0 = 0  # Moment about y
n_0 = 0  # Moment about z

# Parameters
mass = 1          # Mass in kg
J_x = 1        # Moment of Inertia about x in kg*m^2
J_y = 1        # Moment of Inertia about y in kg*m^2
J_z = 2        # Moment of Inertia about z in kg*m^2
J_xy = 0       # Product of Inertia in kg*m^2 (Assumed to be 0 for symmetric body)
J_yz = 0       # Product of Inertia in kg*m^2 (Assumed to be 0 for symmetric body)
J_xz = 0       # Product of Inertia in kg*m^2

# Initial Conditions
## Translational Kinematics State
p_n_0 = 0
p_e_0 = 0
p_d_0 = 100
## Rotational Kinematics State
phi_0 = 0
theta_0 = 0
psi_0 = 0
## Translational Dynamics State
u_0 = 10
v_0 = 0
w_0 = 0
## Rotational Dynamics State
p_0 = 0
q_0 = 0
r_0 = 0

x_IC = np.array([p_n_0, p_e_0, p_d_0, phi_0, theta_0, psi_0, u_0, v_0, w_0, p_0, q_0, r_0]) # Initial Conditions Vector
u_IC = np.array([f_x_0, f_y_0, f_z_0, l_0, m_0, n_0])

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

def f(t, x, u_i):
    p_n, p_e, p_d, phi, theta, psi, u, v, w, p, q, r = x
    f_x, f_y, f_z, l, m, n = u_i
    
    
    # Rotation Matrices
    pos_rot_matrix = np.array([[np.cos(theta)*np.cos(psi), np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi)],
                               [np.cos(theta)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.cos(phi)*np.sin(theta)*np.sin(psi) - np.sin(phi)*np.cos(psi)],
                               [-np.sin(theta), np.sin(phi)*np.cos(theta), np.cos(phi)*np.cos(theta)]])
    euler_rot_matrix = np.array([[1, np.sin(phi)*np.tan(theta), np.cos(phi)*np.tan(theta)],
                                [0, np.cos(phi), -np.sin(phi)],
                                [0, np.sin(phi)/np.cos(theta), np.cos(phi)/np.cos(theta)]])

    pos_dot = np.dot(pos_rot_matrix, np.array([u, v, w])) # From Slide 15
    V_dot = np.array([r*v-q*w, p*w-r*u, q*u-p*v]) + (1/mass)*np.array([f_x, f_y, f_z]) # From Slide 15
    euler_dot = np.dot(euler_rot_matrix, np.array([p, q, r])) # From Slide 15
    W_dot = np.array([Gamma1*p*q-Gamma2*q*r, Gamma5*p*r-Gamma6*(p**2-r**2), Gamma7*p*q-Gamma1*q*r]) + np.array([Gamma3*l+Gamma4*n, (1/J_y)*m, Gamma4*l+Gamma8*n]) # From Slide 15

    return np.concatenate((pos_dot, V_dot, euler_dot, W_dot))

# Time Stepping Parameters
dt = 0.1; n = 100

t = 0; x = x_IC; u_i = u_IC
integrator = intg.RK4(dt, f)
t_history = [0]
x_history = [x]
for i in range(n):
    
    x = integrator.step(t, x, u_i)
    t = (i+1) * dt

    t_history.append(t)
    x_history.append(x)

p_n_values = list(map(lambda x: x[0], x_history))
p_e_values = list(map(lambda x: x[1], x_history))
p_d_values = list(map(lambda x: x[2], x_history))
phi_values = list(map(lambda x: x[3], x_history))
theta_values = list(map(lambda x: x[4], x_history))
psi_values = list(map(lambda x: x[5], x_history))
u_values = list(map(lambda x: x[6], x_history))
v_values = list(map(lambda x: x[7], x_history))
w_values = list(map(lambda x: x[8], x_history))
p_values = list(map(lambda x: x[9], x_history))
q_values = list(map(lambda x: x[10], x_history))
r_values = list(map(lambda x: x[11], x_history))

intg.__doc__
plt.figure()
plt.title('pn')
plt.plot(t_history, p_n_values)
#plt.legend(["pn", "pe", "pd"])
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(p_n_values, p_e_values, p_d_values)
ax.set_xlabel('North')
ax.set_ylabel('East')
ax.set_zlabel('Down')
plt.show()