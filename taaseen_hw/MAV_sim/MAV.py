# In[ ]:
# Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

import integrators as intg
from parameters import *
from inputs import *

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

x_IC = np.array([p_n_0, p_e_0, p_d_0, u_0, v_0, w_0, phi_0, theta_0, psi_0, p_0, q_0, r_0]) # Initial Conditions Vector
u_IC = np.array([u_w0, v_w0, w_w0])

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
    u_w, v_w, w_w = u_i
    
    # State Variable Matrices
    M_ned = np.array([p_n, p_e, p_d])
    M_vg = np.array([u, v, w])
    M_euler = np.array([phi, theta, psi])
    M_angular = np.array([p, q, r])
    
    # Body to Aero Frame
    u_r = u - u_w
    v_r = v - v_w
    w_r = w - w_w
    V_ba = np.array([u_r, v_r, w_r])
    V_a = np.sqrt(u_r**2 + v_r**2 + w_r**2)
    
    # Aerodynamic Calculations
    alpha = np.arctan(w_r/u_r)
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
    J = 2*np.pi*V_a/(Omega*D_prop)
    C_T = C_T2*J**2+C_T1*J+C_T0
    C_Q = C_Q2*J**2+C_Q1*J+C_Q0
    T_p = rho*(Omega/(2*np.pi))**2*D_prop**4*C_T
    Q_p = rho*(Omega/(2*np.pi))**2*D_prop**5*C_Q
    
    # Control Surface Variables
    delta_a = p
    delta_e = -q
    delta_r = -r
    
    # Force Matrix
    force = (np.array([-mass*g*np.sin(theta), 
                       mass*g*np.cos(theta)*np.sin(phi), 
                       mass*g*np.cos(theta)*np.cos(phi)])
       + ((1/2)*rho*V_a**2*S)*np.array([C_X+C_X_q*c/(2*V_a)*q, 
                                        C_Y_0+C_Y_beta*beta+C_Y_p*b/(2*V_a)*p+C_Y_r*b/(2*V_a)*r, 
                                        C_Z+C_Z_q*c/(2*V_a)*q])
       + ((1/2)*rho*V_a**2*S)*np.array([C_X_delta_e*delta_e, 
                                        C_Y_delta_a*delta_a+C_Y_delta_r*delta_r, 
                                        C_Z_delta_e*delta_e])
       + np.array([T_p, 
                   0, 
                   0]))

    # Moment Matrix
    moment = ((1/2)*rho*V_a**2*S*np.array([b*(C_l_0+C_l_beta*beta+C_l_p*b/(2*V_a)*p+C_l_r*b/(2*V_a)*r), 
                                           c*(C_m_0+C_m_alpha*alpha+C_m_q*c/(2*V_a)*q),
                                           b*(C_n_0+C_n_beta*beta+C_n_p*b/(2*V_a)*p+C_n_r*b/(2*V_a)*r)])
            + (1/2)*rho*V_a**2*S*np.array([b*(C_l_delta_a*delta_a+C_l_delta_r*delta_r),
                                           c*(C_m_delta_e*delta_e),
                                           b*(C_n_delta_a*delta_a+C_n_delta_r*delta_r)])
            + np.array([Q_p, 
                        0, 
                        0]))
    
    # Rotation Matrices
    pos_rot_matrix = np.array([[np.cos(theta)*np.cos(psi), np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi)],
                               [np.cos(theta)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.cos(phi)*np.sin(theta)*np.sin(psi) - np.sin(phi)*np.cos(psi)],
                               [-np.sin(theta), np.sin(phi)*np.cos(theta), np.cos(phi)*np.cos(theta)]])
    euler_rot_matrix = np.array([[1, np.sin(phi)*np.tan(theta), np.cos(phi)*np.tan(theta)],
                                [0, np.cos(phi), -np.sin(phi)],
                                [0, np.sin(phi)/np.cos(theta), np.cos(phi)/np.cos(theta)]])

    pos_dot = np.dot(pos_rot_matrix, M_vg)
    V_dot = np.array([r*v-q*w, p*w-r*u, q*u-p*v]) + (1/mass)*force
    euler_dot = np.dot(euler_rot_matrix, M_angular)
    W_dot = np.linalg.inv(J_i) @ (np.cross(-M_angular, (J_i @ M_angular)) + moment)
    
    return np.concatenate((pos_dot, V_dot, euler_dot, W_dot))

# Time Stepping Parameters
dt = 0.01; num = 10

t = 0; x = x_IC; u_i = u_IC
integrator = intg.RK4(dt, f)
t_history = [0]
x_history = [x]
for i in range(num):
    
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