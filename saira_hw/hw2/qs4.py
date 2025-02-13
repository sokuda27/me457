import numpy as np

# Inputs
# forces applied to drone
fx=1 #force in x
fy=1
fz=1
F_b = np.array([[fx], [fy], [fz]])

#moments applied to drone
l=1 #moment about x 
m=1
n=1
M_b = np.array([[l],[m],[n]])

# Parameters
mass=1
#moments of inertia
Jx=1; Jy=1; Jz=1
Jxy=0; Jyz=0; Jxz=0
J=np.array([[Jx, -Jxy, -Jxz],
   [-Jxy, Jy, -Jyz],
   [-Jxz, -Jyz, Jz]])

# Initial Conditions
u=1; v=1; w=1
v_linear = np.array([[u],[v],[w]])
p=1; q=1; r=1
w_angular = np.array([[p],[q],[r]])
w_bb = np.array([[0, -r, q],
        [r, 0, -p],
        [-q, p, 0]])

v_linear_sol = (1/mass**2)*np.dot(-w_bb, v_linear)+F_b
print('udot, vdot, wdot:')
print(v_linear_sol)
w_angular_sol = np.dot(np.linalg.inv(J), (M_b-np.dot(w_bb,np.dot(J,w_angular))))
print('pdot, qdot, rdot')
print(w_angular_sol)
