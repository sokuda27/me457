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

# Mass and Inertia Parameters
mass = 1       # Mass in kg
g = 9.81       # Gravitational Acceleration in m/s^2
rho = 1.268    # Air Density in kg/m^3
J_x = 1        # Moment of Inertia about x in kg*m^2
J_y = 1        # Moment of Inertia about y in kg*m^2
J_z = 2        # Moment of Inertia about z in kg*m^2
J_xy = 0       # Product of Inertia in kg*m^2 (Assumed to be 0 for symmetric body)
J_yz = 0       # Product of Inertia in kg*m^2 (Assumed to be 0 for symmetric body)
J_xz = 0       # Product of Inertia in kg*m^2

# Aerodynamic Parameters
C_D_0 = 0.043
C_L_0 = 0.23
C_D_alpha = 0.030
C_L_alpha = 5.61
C_D_q = 0
C_L_q = 7.95
C_D_delta_e = 0.0135
C_L_delta_e = 0.13
C_Y_0 = 0
C_Y_delta_a = 0.075
C_Y_delta_r = 0.19
C_Y_beta = -0.83
C_Y_p = 0
C_Y_r = 0

# Physical Parameters
S = 0.55 # Wing Area in m^2
b = 2.9  # Wing Span in m
c = 0.19 # Mean Aerodynamic Chord in m
D = 0.508 # Diameter of Propeller in m
V_max = 44.4 # Maximum Velocity in m/s
K_V = 0.0659 # V*s/rad
Omega = V_max/K_V # Propellor Speed in rad/s
C_T0 = 0.09357
C_T1 = -0.06044
C_T2 = -0.1079
C_Q0 = 0.005230
C_Q1 = 0.004970
C_Q2 = -0.01664

# Wind Parameters
u_w = 0
v_w = 0
w_w = 0

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

x_IC = np.array([p_n_0, p_e_0, p_d_0, u_0, v_0, w_0, phi_0, theta_0, psi_0, p_0, q_0, r_0]) # Initial Conditions Vector
u_IC = np.array([l_0, m_0, n_0])

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
    p_n, p_e, p_d, u, v, w, phi, theta, psi, p, q, r = x
    l, m, n = u_i
    
    # Body to Aero Frame
    u_r = u - u_w
    v_r = v - v_w
    w_r = w - w_w
    V_ba = np.array([u_r, v_r, w_r])
    V_a = np.sqrt(u_r**2 + v_r**2 + w_r**2)
    
    # Aerodynamic Calculations
    alpha = np.arctan2(w_r, u_r)
    beta = np.arcsin(v_r/(np.sqrt(u_r**2 + v_r**2 + w_r**2)))
    C_D = C_D_0 + C_D_alpha*alpha
    C_L = C_L_0 + C_L_alpha*alpha
    C_X = -C_D*np.cos(alpha) + C_L*np.sin(alpha)
    C_X_q = -C_D_q*np.cos(alpha) + C_L_q*np.sin(alpha)
    C_X_delta_e = -C_D_delta_e*np.cos(alpha) + C_L_delta_e*np.sin(alpha)
    C_Z = -C_D*np.sin(alpha) - C_L*np.cos(alpha)
    C_Z_q = -C_D_q*np.sin(alpha) - C_L_q*np.cos(alpha)
    C_Z_delta_e = -C_D_delta_e*np.sin(alpha) - C_L_delta_e*np.cos(alpha)

    # Torque and Thrust Calculations
    J = 2*np.pi*V_a/(Omega*D)
    C_T = C_T2*J**2+C_T1*J+C_T0
    C_Q = C_Q2*J**2+C_Q1*J+C_Q0
    T_p = rho*(Omega/(2*np.pi))**2*D**4*C_T
    
    # Force Matrix
    f = (np.array([-mass*g*np.sin(theta), mass*g*np.cos(theta)*np.sin(phi), mass*g*np.cos(theta)*np.cos(phi)])
       + ((1/2)*rho*V_a**2*S)*np.array([C_X+C_X_q*c/(2*V_a)*q, C_Y_0+C_Y_beta*beta+C_Y_p*b/(2*V_a*p+C_Y_r*b/(2*V_a)*r), C_Z+C_Z_q*c/(2*V_a)*q])
       + ((1/2)*rho*V_a**2*S)*np.array([C_X_delta_e*(-q), C_Y_delta_a*p+C_Y_delta_r*(-r), C_Z_delta_e*(-q)])  
       + np.array([T_p, 0, 0]))
    
    # Rotation Matrices
    pos_rot_matrix = np.array([[np.cos(theta)*np.cos(psi), np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi)],
                               [np.cos(theta)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.cos(phi)*np.sin(theta)*np.sin(psi) - np.sin(phi)*np.cos(psi)],
                               [-np.sin(theta), np.sin(phi)*np.cos(theta), np.cos(phi)*np.cos(theta)]])
    euler_rot_matrix = np.array([[1, np.sin(phi)*np.tan(theta), np.cos(phi)*np.tan(theta)],
                                [0, np.cos(phi), -np.sin(phi)],
                                [0, np.sin(phi)/np.cos(theta), np.cos(phi)/np.cos(theta)]])

    pos_dot = np.dot(pos_rot_matrix, np.array([u, v, w])) # From Slide 15
    V_dot = np.array([r*v-q*w, p*w-r*u, q*u-p*v]) + (1/mass)*f # From Slide 15
    euler_dot = np.dot(euler_rot_matrix, np.array([p, q, r])) # From Slide 15
    W_dot = np.array([Gamma1*p*q-Gamma2*q*r, Gamma5*p*r-Gamma6*(p**2-r**2), Gamma7*p*q-Gamma1*q*r]) + np.array([Gamma3*l+Gamma4*n, (1/J_y)*m, Gamma4*l+Gamma8*n]) # From Slide 15

    print(V_a)
    
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

p_n_values = [arr[0] for arr in x_history]
p_e_values = [arr[1] for arr in x_history]
p_d_values = [arr[2] for arr in x_history]
u_values = [arr[3] for arr in x_history]
v_values = [arr[4] for arr in x_history]
w_values = [arr[5] for arr in x_history]
phi_values = [arr[6] for arr in x_history]
theta_values = [arr[7] for arr in x_history]
psi_values = [arr[8] for arr in x_history]
p_values = [arr[9] for arr in x_history]
q_values = [arr[10] for arr in x_history]
r_values = [arr[11] for arr in x_history]


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