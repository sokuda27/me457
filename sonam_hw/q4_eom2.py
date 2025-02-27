import numpy as np
import q4_parameters as param
import integrators_sonam as intg
import matplotlib.pyplot as plt

# Defining pn, pe, pd
def state(t, x, u):
    rot_mat_1 = np.array([[np.cos(x[7])*np.cos(x[8]), np.sin(x[6])*np.sin(x[7])*np.cos(x[8])-np.cos(6)*np.sin(x[8]), np.cos(x[6])*np.sin(x[7])*np.cos(x[8])+np.sin(x[6])*np.sin(x[8])],
                          [np.cos(x[7])*np.cos(x[8]), np.sin(x[6])*np.sin(x[7])*np.sin(x[8])+np.cos(6)*np.cos(x[8]), np.cos(x[6])*np.sin(x[7])*np.sin(x[8])-np.sin(x[6])*np.cos(x[8])],
                          [-np.sin(x[7]), np.sin(x[6])*np.cos(x[7]), np.cos(x[6])*np.cos(x[7])]
                          ])
    pned_dot = rot_mat_1 @ np.array([x[3], x[4], x[5]]).T
    uvw_dot = np.array([[x[11]*x[4]-x[10]*x[5]],
                        [x[9]*x[5]-x[11]*x[3]],
                        [x[10]*x[3]-x[9]*x[4]]]) + (1/param.mass)*np.array([u[0], u[1], u[2]])
    euler_dot = np.array([[1, np.sin(x[6])*np.tan(x[7]), np.cos(x[6])*np.tan(x[7])],
                          [0, np.cos(x[6]), -np.sin(x[6])],
                          [0, np.sin(x[6])/np.cos(x[7]), np.cos(x[6])*np.cos(x[7])]]) @ np.array([x[9], x[10], x[11]]).T
    pqr_dot1 = np.array([[param.j1*x[9]*x[10]-param.j2*x[10]*x[11]],
                        [param.j5*x[9]*x[11]-param.j6*(x[9]**2-x[11]**2)],
                        [param.j7*x[9]*x[10]-param.j1*x[10]*x[11]]]) + param.pqr_1
    return np.array([pned_dot[0], pned_dot[1], pned_dot[2], uvw_dot[0], uvw_dot[1], uvw_dot[2], euler_dot[0], euler_dot[1], euler_dot[2], pqr_dot1[0], pqr_dot1[1], pqr_dot1[2]])
    
# inital conditions x
p_n0 = 0; p_e0 = 0; p_d0 = 100
u_0 = 10; v_0 = 0; w_0 = 0
phi_0 = 0; theta_0 = 0; psi_0 = 0
p_0 = 0; q_0 = 0; r_0 = 0

# initial conditions u
Fx_0 = 1; Fy_0 = 0; Fz_0 = 9.8*param.mass
lx_0 = 0; my_0 = 0; nz_0 = 0

# setting up plotting parameters
t = 0
x = np.array([p_n0, p_e0, p_d0, u_0, v_0, w_0, phi_0, theta_0, psi_0, p_0, q_0, r_0])
u_i = np.array([Fx_0, Fy_0, Fz_0, lx_0, my_0, nz_0])
dt = .1; n = 100

integrator = intg.RungeKutta4(dt, state)

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

