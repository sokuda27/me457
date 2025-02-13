import numpy as np
import integrators_sonam as intg
import matplotlib.pyplot as plt

# Inputs
# Forces
fx  = 1
fy = 1
fz = 1
F_b = np.array([[fx], [fy], [fz]])
# Moments
m1 = 1
m2 = 1
m3 = 1
M_b = np.array([[m1], [m2], [m3]])

# Parameters
mass = 11 #kg
Jx = 0.824 #kg-m^2
Jy = 1.135
Jz = 1.759 
Jxz = 0.120
J_bb = [[Jx, 0, Jxz], [0, Jy, 0], [0, 0, Jz]]

# Initial conditions
u = 0
v = 0
w = 0
u = 1
v = 1
w = 1
p = 1
q = 1
r = 1

v_b = [u, v, w].T
w_b = [p, q, r].T
w_bb = np.array([[0, -r, q], [r, 0, -p], [-q, p, 0]])

# Formulas
v_sol = (1/mass**2)

# # Setting up equations
# def omega_p(t, w, u):
#     w_temp = w[0:2]
#     w_dot_p = np.linalg.inv(J_bb)@(M_b - np.cross(w_temp, np.cross(J_bb, w_temp)))
#     return np.array([w[3], w_dot_p[0]])

# # Setting up integrator

# # t = 0; x = np.array([0, 0, 0, 1, 1, 1]); u = 0
# dt = 0.1; n = 50
# integrator = intg.RungeKutta4(dt, omega_p)

# t_history = [0] 
# x_history = [x]
# for i in range(n):
    
#     x = integrator.step(t, x, u)
#     t = (i+1) * dt

#     t_history.append(t)
#     x_history.append(x)

# plt.figure()
# plt.plot(t_history, x_history)
# plt.legend(["theta","omega"])
# plt.show()