import numpy as np

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
J_bb = [[Jx, 0, 0], [0, Jy, 0], [Jxz, 0, 0]]

# Initial conditions
x = 0
y = 0
z = 0
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