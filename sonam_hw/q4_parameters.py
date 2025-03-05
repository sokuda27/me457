import numpy as np

# Inputs
# # Forces
# fx  = 0
# fy = 0
# fz = 0
# F_b = np.array([[fx], [fy], [fz]])
# # Moments
m = 0,; n = 0; l = 0

# Parameters
# mass = 11 #kg
# Jx = 0.824 #kg-m^2
# Jy = 1.135
# Jz = 1.759 
# Jxz = 0.120
# J_bb = [[Jx, 0, Jxz], [0, Jy, 0], [0, 0, Jz]]

mass = 11 #kg
g = 9.81
Jx = 0.824 #kg-m^2
Jy = 1.135
Jz = 1.759
Jxy=0; Jyz=0; Jxz=0
Jyx=0; Jzy=0; Jzx=0
J_bb = [[Jx, Jxy, Jxz], [Jyx, Jy, Jyz], [Jzx, Jzy, Jz]]

jj = Jx*Jz-Jxz**2
j1 = (Jxz*(Jx-Jy+Jz))/jj
j2 = (Jz*(Jz-Jy)+Jxz**2)/jj
j3 = Jz/jj
j4 = Jxz/jj
j5 = (Jz - Jx)/Jy
j6 = Jxz/Jy
j7 = ((Jx - Jy)*Jx + Jxz**2)/jj
j8 = Jx/jj

# pqr_1 = np.array([[j3*l+j4*n],
#                   [0], # m/Jy ddnt work?
#                   [j4*l+j8*n]])

rho = 1.268 #kg/m^3
S = 0.55 #m^2
b = 2.9 #m
c = 0.19
v_max = 44.4 #V
kv = 0.0659 #Vs/rad
omg = v_max/kv
D = 0.508 #m

cl_0 = 0.23
cl_alp = 5.61
cd_0 = 0.043
cd_alp = 0.030

cd_q = 0
cl_q = 7.95

cd_de = 0.0135
cl_de = 0.13

ct2 = -0.1079
ct1 = -0.06044
ct0 = 0.09357

cq1 = 0.004970
cq2 = -0.01664
cq0 = 0.005230

cy_0 = 0 
cy_be = -0.83
cy_p = 0; cy_r = 0
cy_da = 0.075
cy_dr = 0.19

cl_p = -0.51
cl_be = -0.13
cl_r = 0.25
cl_da = 0.17
cl_dr = 0.0024

cm_0 = 0.0135
cm_alp = -2.74
cm_q = -38.21
cm_de = -0.99

cn_0 = 0
cn_p = -0.069
cn_be = 0.073
cn_r = -0.095
cn_da = -0.011
cn_dr = -0.069

def aero_coeff(alpha, va):
    cd = cl_0 + cl_alp*alpha
    cl = cd_0 + cd_alp*alpha
    cx = -cd*np.cos(alpha) + cl*np.sin(alpha)
    cx_q = -cd_q*np.cos(alpha) + cl_q*np.sin(alpha)
    cx_de = -cd_de*np.cos(alpha) + cl_de*np.sin(alpha)
    cz = -cd*np.sin(alpha) - cl*np.cos(alpha)
    cz_q = -cd_q*np.sin(alpha) - cl_q*np.cos(alpha)
    cz_de = -cd_de*np.sin(alpha) - cl_de*np.cos(alpha)

    j_aero = 2*np.pi*va/(omg*D)
    ct = ct2*j_aero**2 + ct1*j_aero + ct0
    cq = cq2*j_aero**2 + ct1*j_aero + ct0
    tp = rho*(omg/(2*np.pi))**2*D**4*ct
    qp = rho*(omg/(2*np.pi))**2*D**5*cq

    return np.array([cx, cx_q, cx_de, cz, cz_q, cz_de, tp, qp])

def f(t, x, u):
    u_r = x[3] - u[0]
    v_r = x[4] - u[1]
    w_r = x[5] - u[2]
    alp = np.arctan(w_r/u_r)
    v_a = np.sqrt(u_r**2 + v_r**2 + w_r**2)
    beta = np.arcsin(v_r/v_a)
    coeffs = aero_coeff(alp, v_a)

    fx = -mass*g*np.sin(x[7]) + 0.5*v_a**2*S*(coeffs[0] + coeffs[1]*c*x[10]/(2*v_a) - coeffs[2]*x[10] + coeffs[6])
    fy = mass*g*np.cos(x[7])*np.sin(x[6]) + 0.5*v_a**2*S*(cy_0 + cy_be*beta + cy_p*b/(2*v_a)*x[9] + cy_r*b/(2*v_a)*x[11] + cy_da*x[9] - cy_dr*x[11])
    fz = mass*g*np.cos(x[7])*np.cos(x[6]) + 0.5*v_a**2*S*(coeffs[3] + coeffs[4]*c/(2*v_a)*x[10] - coeffs[5]*x[10])

    l = 0.5*v_a**2*S*b*(cl_0 + cl_be*beta + cl_p*(b/(2*v_a))*x[9] + cl_r*(b/(2*v_a))*x[11] + cl_da*x[9] + cl_dr*(-x[11])) + coeffs[7]
    m = 0.5*v_a**2*S*c*(cm_0 + cm_alp*alp + cm_q*(c/(2*v_a))*x[10] + cm_de*x[10])
    n = 0.5*v_a**2*S*b*(cn_0 + cn_be*beta + cn_p*(b/(2*v_a))*x[9] + cn_r*(b/(2*v_a))*x[11] + cn_da*x[9] + cn_dr*(-x[11]))

    pn_dot = np.cos(x[7])*np.cos(x[8])*x[3] + (np.sin(x[6])*np.sin(x[7])*np.cos(x[8])-np.cos(6)*np.sin(x[8]))*x[4] + (np.cos(x[6])*np.sin(x[7])*np.cos(x[8])+np.sin(x[6])*np.sin(x[8]))*x[5]
    pe_dot = np.cos(x[7])*np.sin(x[8])*x[3] + (np.sin(x[6])*np.sin(x[7])*np.sin(x[8])+np.cos(6)*np.cos(x[8]))*x[4] + (np.cos(x[6])*np.sin(x[7])*np.sin(x[8])-np.sin(x[6])*np.cos(x[8]))*x[5]
    pd_dot = (-np.sin(x[7]))*x[3] + np.sin(x[6])*np.cos(x[7])*x[4] + np.cos(x[6])*np.cos(x[7])*x[5]
    u_dot = x[11]*x[4]-x[10]*x[5] + (1/mass)*fx
    v_dot = x[9]*x[5]-x[11]*x[3] + (1/mass)*fy
    w_dot = x[10]*x[3]-x[9]*x[4] + (1/mass)*fz
    phi_dot = x[9] + np.sin(x[6])*np.tan(x[7])*x[10] + np.cos(x[6])*np.tan(x[7])*x[11]
    theta_dot = np.cos(x[6])*x[10] - np.sin(x[6])*x[11]
    psi_dot = np.sin(x[6])/np.cos(x[7])*x[10] + np.cos(x[6])*np.cos(x[7])*x[11]
    p_dot = j1*x[9]*x[10]-j2*x[10]*x[11] + j3*l+j4*n
    q_dot = j5*x[9]*x[11]-j6*(x[9]**2 - x[11]**2) + m/Jy
    r_dot = j7*x[9]*x[10]-j1*x[10]*x[11] + j4*l+j8*n
    return np.array([pn_dot, pe_dot, pd_dot, u_dot, v_dot, w_dot, phi_dot, theta_dot, psi_dot, p_dot, q_dot, r_dot])
