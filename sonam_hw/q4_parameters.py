import numpy as np

# Inputs
# Forces
fx  = 0
fy = 0
fz = 0
F_b = np.array([[fx], [fy], [fz]])
# Moments
m1 = 1
m2 = 1
m3 = 1
M_b = np.array([[m1], [m2], [m3]])

# Parameters
# mass = 11 #kg
# Jx = 0.824 #kg-m^2
# Jy = 1.135
# Jz = 1.759 
# Jxz = 0.120
# J_bb = [[Jx, 0, Jxz], [0, Jy, 0], [0, 0, Jz]]

mass = 1 #kg
Jx = 1 #kg-m^2
Jy = 1
Jz = 2
Jxy=0; Jyz=0; Jxz=0
Jyx=0; Jzy=0; Jzx=0
J_bb = [[Jx, Jxy, Jxz], [Jyx, Jy, Jyz], [Jzx, Jzy, Jz]]

jj = Jx*Jz-Jxz^2
j1 = (Jxz*(Jx-Jy+Jz))/jj
j2 = (Jz*(Jz-Jy)+Jxz^2)/jj
j3 = Jz/jj
j4 = Jxz/jj
j5 = (Jz - Jx)/Jy
j6 = Jxz/Jy
j7 = ((Jx - Jy)*Jx + Jxz^2)/jj
j8 = Jx/jj

m = 0, n = 0, l = 0
pqr_1 = np.array([[j3*l+j4*n],
                  [(1/Jy)*m],
                  [j4*l+j8*n]])
