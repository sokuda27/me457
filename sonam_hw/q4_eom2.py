import numpy as np
import q4_parameters as param

# Defining pn, pe, pd
def state(t, x, u):
    rot_mat_1 = np.array([[np.cos(x[7])*np.cos(x[8]), np.sin(x[6])*np.sin(x[7])*np.cos(x[8])-np.cos(6)*np.sin(x[8]), np.cos(x[6])*np.sin(x[7])*np.cos(x[8])+np.sin(x[6])*np.sin(x[8])],
                          [np.cos(x[7])*np.cos(x[8]), np.sin(x[6])*np.sin(x[7])*np.sin(x[8])+np.cos(6)*np.cos(x[8]), np.cos(x[6])*np.sin(x[7])*np.sin(x[8])-np.sin(x[6])*np.cos(x[8])],
                          [-np.sin(x[7]), np.sin(x[6])*np.cos(x[7]), np.cos(x[6])*np.cos(x[7])]
                          ])
    pned_dot = rot_mat_1 @ x[3:6]
    uvw_dot = np.array([[x[12]*x[4]-x[10]*x[5]],
                        [x[9]*x[5]-x[11]*x[3]],
                        [x[10]*x[3]-x[9]*x[4]]]) + (1/param.mass)*(u[0], u[1], u[2])
    euler_dot = np.array([[1, np.sin(x[6])*np.tan(x[7]), np.cos(x[6])*np.tan(x[7])],
                          [0, np.cos(x[6]), -np.sin(x[6])],
                          [0, np.sin(x[6])/np.cos(x[7]), np.cos(x[6])*np.cos(x[7])]]) @ x[9:]
    pqr_dot1 = np.array([[param.j1*x[9]*x[10]-param.j2*x[10]*x[11]],
                        [param.j5*x[9]*x[11]-param.j6*(x[9]**2-x[11]**2)],
                        [param.j7*x[9]*x[10]-param.j1*x[10]*x[11]]]) + param.pqr_1