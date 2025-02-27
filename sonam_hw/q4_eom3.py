import numpy as np
import q4_parameters as param
import integrators_sonam as intg
import matplotlib.pyplot as plt

# Defining state 
def f(t, x, u):
    u_r = x[3] - u[0]
    v_r = x[4] - u[1]
    w_r = x[5] - u[2]
    alp = np.arctan(w_r/u_r)
    v_a = np.sqrt(u_r**2 + v_r**2 + w_r**2)
    beta = np.arcsin(v_r/v_a)
    coeffs = param.aero_coeff(alp, v_a)

    fx = -param.mass*param.g*np.sin(x[7]) + 0.5*v_a**2*param.S*(coeffs[0] + coeffs[1]*param.c*x[10]/(2*v_a) - coeffs[2]*x[10] + coeffs[6])
    fy = param.mass*param.g*np.cos(x[7])*np.sin(x[6]) + 0.5*v_a**2*param.S*(param.cy_0 + param.cy_be*beta + param.cy_p*param.b/(2*v_a)*x[9] + param.cy_r*param.b/(2*v_a)*x[11] + param.cy_da*x[9] - param.cy_dr*x[11])
    fz = param.mass*param.g*np.cos(x[7])*np.cos(x[6]) + 0.5*v_a**2*param.S*(coeffs[3] + coeffs[4]*param.c/(2*v_a)*x[10] - coeffs[5]*x[10])

    l = 0.5*v_a**2*param.S*param.b*(param.cl_0 + param.cl_be*beta + param.cl_p*(param.b/(2*v_a))*x[9] + param.cl_r*(param.b/(2*v_a))*x[11] + param.cl_da*x[9] + param.cl_dr*(-x[11])) + coeffs[7]
    m = 0.5*v_a**2*param.S*param.c*(param.cm_0 + param.cm_alp*alp + param.cm_q*(param.c/(2*v_a))*x[10] + param.cm_de*x[10])
    n = 0.5*v_a**2*param.S*param.b*(param.cn_0 + param.cn_be*beta + param.cn_p*(param.b/(2*v_a))*x[9] + param.cn_r*(param.b/(2*v_a))*x[11] + param.cn_da*x[9] + param.cn_dr*(-x[11]))

    pn_dot = np.cos(x[7])*np.cos(x[8])*x[3] + (np.sin(x[6])*np.sin(x[7])*np.cos(x[8])-np.cos(6)*np.sin(x[8]))*x[4] + (np.cos(x[6])*np.sin(x[7])*np.cos(x[8])+np.sin(x[6])*np.sin(x[8]))*x[5]
    pe_dot = np.cos(x[7])*np.sin(x[8])*x[3] + (np.sin(x[6])*np.sin(x[7])*np.sin(x[8])+np.cos(6)*np.cos(x[8]))*x[4] + (np.cos(x[6])*np.sin(x[7])*np.sin(x[8])-np.sin(x[6])*np.cos(x[8]))*x[5]
    pd_dot = (-np.sin(x[7]))*x[3] + np.sin(x[6])*np.cos(x[7])*x[4] + np.cos(x[6])*np.cos(x[7])*x[5]
    u_dot = x[11]*x[4]-x[10]*x[5] + (1/param.mass)*fx
    v_dot = x[9]*x[5]-x[11]*x[3] + (1/param.mass)*fy
    w_dot = x[10]*x[3]-x[9]*x[4] + (1/param.mass)*fz
    phi_dot = x[9] + np.sin(x[6])*np.tan(x[7])*x[10] + np.cos(x[6])*np.tan(x[7])*x[11]
    theta_dot = np.cos(x[6])*x[10] - np.sin(x[6])*x[11]
    psi_dot = np.sin(x[6])/np.cos(x[7])*x[10] + np.cos(x[6])*np.cos(x[7])*x[11]
    p_dot = param.j1*x[9]*x[10]-param.j2*x[10]*x[11] + param.j3*l+param.j4*n
    q_dot = param.j5*x[9]*x[11]-param.j6*(x[9]**2 - x[11]**2) + m/param.Jy
    r_dot = param.j7*x[9]*x[10]-param.j1*x[10]*x[11] + param.j4*l+param.j8*n
    return np.array([pn_dot, pe_dot, pd_dot, u_dot, v_dot, w_dot, phi_dot, theta_dot, psi_dot, p_dot, q_dot, r_dot])

# inital conditions x
p_n0 = 1; p_e0 = 0; p_d0 = 2
u_0 = 1; v_0 = 1; w_0 = 0
phi_0 = 0; theta_0 = 0; psi_0 = 0
p_0 = 0; q_0 = 0; r_0 = 0

x = np.array([p_n0, p_e0, p_d0, u_0, v_0, w_0, phi_0, theta_0, psi_0, p_0, q_0, r_0])

# initial conditions u
# Fx_0 = 1; Fy_0 = 0; Fz_0 = 9.8*param.mass
# lx_0 = 0; my_0 = 0; nz_0 = 0

# wind input
u_w = 10; v_w = 10; w_w = 0

ui = np.array([u_w, v_w, w_w])

# setting up integrator
dt = 0.01; n = 10; t = 0

integrator = intg.RungeKutta4(dt, f)

t_history = [0]
x_history = [x]
for i in range(n):
    
    x = integrator.step(t, x, ui)
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