"""
compute_ss_model
    - Chapter 5 assignment for Beard & McLain, PUP, 2012
    - Update history:  
        2/4/2019 - RWB
"""
import numpy as np
from scipy.optimize import minimize
from tools.rotations import euler_to_quaternion, quaternion_to_euler
import parameters.aerosonde_parameters as MAV
from parameters.simulation_parameters import ts_simulation as Ts
from message_types.msg_delta import MsgDelta


def compute_model(mav, trim_state, trim_input):
    # Note: this function alters the mav private variables
    A_lon, B_lon, A_lat, B_lat = compute_ss_model(mav, trim_state, trim_input)
    Va_trim, alpha_trim, theta_trim, a_phi1, a_phi2, a_theta1, a_theta2, a_theta3, \
    a_V1, a_V2, a_V3 = compute_tf_model(mav, trim_state, trim_input)

    # write transfer function gains to file
    file = open('models/model_coef.py', 'w')
    file.write('import numpy as np\n')
    file.write('x_trim = np.array([[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]]).T\n' %
               (trim_state.item(0), trim_state.item(1), trim_state.item(2), trim_state.item(3),
                trim_state.item(4), trim_state.item(5), trim_state.item(6), trim_state.item(7),
                trim_state.item(8), trim_state.item(9), trim_state.item(10), trim_state.item(11),
                trim_state.item(12)))
    file.write('u_trim = np.array([[%f, %f, %f, %f]]).T\n' %
               (trim_input.elevator, trim_input.aileron, trim_input.rudder, trim_input.throttle))
    file.write('Va_trim = %f\n' % Va_trim)
    file.write('alpha_trim = %f\n' % alpha_trim)
    file.write('theta_trim = %f\n' % theta_trim)
    file.write('a_phi1 = %f\n' % a_phi1)
    file.write('a_phi2 = %f\n' % a_phi2)
    file.write('a_theta1 = %f\n' % a_theta1)
    file.write('a_theta2 = %f\n' % a_theta2)
    file.write('a_theta3 = %f\n' % a_theta3)
    file.write('a_V1 = %f\n' % a_V1)
    file.write('a_V2 = %f\n' % a_V2)
    file.write('a_V3 = %f\n' % a_V3)
    file.write('A_lon = np.array([\n    [%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f]])\n' %
    (A_lon[0][0], A_lon[0][1], A_lon[0][2], A_lon[0][3], A_lon[0][4],
     A_lon[1][0], A_lon[1][1], A_lon[1][2], A_lon[1][3], A_lon[1][4],
     A_lon[2][0], A_lon[2][1], A_lon[2][2], A_lon[2][3], A_lon[2][4],
     A_lon[3][0], A_lon[3][1], A_lon[3][2], A_lon[3][3], A_lon[3][4],
     A_lon[4][0], A_lon[4][1], A_lon[4][2], A_lon[4][3], A_lon[4][4]))
    file.write('B_lon = np.array([\n    [%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f]])\n' %
    (B_lon[0][0], B_lon[0][1],
     B_lon[1][0], B_lon[1][1],
     B_lon[2][0], B_lon[2][1],
     B_lon[3][0], B_lon[3][1],
     B_lon[4][0], B_lon[4][1],))
    file.write('A_lat = np.array([\n    [%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f]])\n' %
    (A_lat[0][0], A_lat[0][1], A_lat[0][2], A_lat[0][3], A_lat[0][4],
     A_lat[1][0], A_lat[1][1], A_lat[1][2], A_lat[1][3], A_lat[1][4],
     A_lat[2][0], A_lat[2][1], A_lat[2][2], A_lat[2][3], A_lat[2][4],
     A_lat[3][0], A_lat[3][1], A_lat[3][2], A_lat[3][3], A_lat[3][4],
     A_lat[4][0], A_lat[4][1], A_lat[4][2], A_lat[4][3], A_lat[4][4]))
    file.write('B_lat = np.array([\n    [%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f]])\n' %
    (B_lat[0][0], B_lat[0][1],
     B_lat[1][0], B_lat[1][1],
     B_lat[2][0], B_lat[2][1],
     B_lat[3][0], B_lat[3][1],
     B_lat[4][0], B_lat[4][1],))
    file.write('Ts = %f\n' % Ts)
    file.close()


def compute_tf_model(mav, trim_state, trim_input):
    # trim values
    mav._state = trim_state
    mav._update_velocity_data()
    Va_trim = mav._Va
    alpha_trim = mav._alpha
    phi, theta_trim, psi = quaternion_to_euler(trim_state[6:10])
    delta_trim = mav._delta

    ###### TODO ######
    # define transfer function constants
    a_phi1 = -0.5 * MAV.rho * Va_trim**2 *MAV.S_wing * MAV.b**2 * MAV.C_p_p /(2*Va_trim)
    a_phi2 = 0.5 * MAV.rho * Va_trim**2 *MAV.S_wing * MAV.b * MAV.C_p_delta_a
    a_theta1 = - (MAV.rho * Va_trim**2 * MAV.c * MAV.S_wing * MAV.C_m_q * MAV.c)/(2 * MAV.Jy * 2 * Va_trim)
    a_theta2 = - (MAV.rho * Va_trim**2 * MAV.c * MAV.S_wing * MAV.C_m_alpha)/(2 * MAV.Jy)
    a_theta3 = - (MAV.rho * Va_trim**2 * MAV.c * MAV.S_wing * MAV.C_m_delta_e)/(2 * MAV.Jy)

    # Compute transfer function coefficients using new propulsion model
    a_V1 = MAV.rho * Va_trim * MAV.S_wing * (MAV.C_D_0 + MAV.C_D_alpha*alpha_trim + MAV.C_D_delta_e*delta_trim) 
    a_V2 = 0
    a_V3 = 0

    return Va_trim, alpha_trim, theta_trim, a_phi1, a_phi2, a_theta1, a_theta2, a_theta3, a_V1, a_V2, a_V3


def compute_ss_model(mav, trim_state, trim_input):
    x_euler = euler_state(trim_state)
    
    ##### TODO #####
    A = df_dx(mav, x_euler, trim_input)
    B = df_du(mav, x_euler, trim_input)
    
        # States Assignments
    p_n = 0
    p_e = 1
    h = 2
    u = 3
    v = 4
    w = 5
    phi = 6 
    theta = 7
    psi = 8
    p = 9
    q = 10
    r = 11
    
    # Inputs Assignments
    delta_e = 0
    delta_a = 1
    delta_r = 2
    delta_t = 3
    
    # extract longitudinal states (u, w, q, theta, pd)
    A_lon = np.zeros((5,5))
    A_lon = np.array([[A[u,u], A[u,w], A[u,q], A[u,theta], A[u,h]],
                      [A[w,u], A[w,w], A[w,q], A[w,theta], A[w,h]],
                      [A[q,u], A[q,w], A[q,q], A[q,theta], A[q,h]],
                      [A[theta,u], A[theta,w], A[theta,q], A[theta,theta], A[theta,h]],
                      [A[h,u], A[h,w], A[h,q], A[h,theta], A[h,h]]])
    B_lon = np.zeros((5,2))
    B_lon = np.array([[[u,delta_e], [u,delta_t]],
                      [[w, delta_e], [w,delta_t]],
                      [[q,delta_e], [q,delta_t]],
                      [[theta,delta_e], [theta,delta_t]],
                      [[h,delta_e], [h,delta_t]]])

    # extract lateral states (v, p, r, phi, psi)
    # A_lat = np.zeros((5,5))
    A_lat = np.array([[A[v,v], A[v,p], A[v,r], A[v,phi], A[v,psi]],
                      [A[p,v], A[p,p], A[p,r], A[p,phi], A[p,psi]],
                      [A[r,v], A[r,p], A[r,r], A[r,phi], A[r,psi]],
                      [A[phi,v], A[phi,p], A[phi,r], A[phi,phi], A[phi,psi]],
                      [A[psi,v], A[psi,p], A[psi,r], A[psi,phi], A[psi,psi]]])
    # B_lat = np.zeros((5,2))
    B_lat = np.array([[B[v,delta_a], B[v,delta_r]],
                      [B[p, delta_a], B[p, delta_r]],
                      [B[r, delta_a], B[r, delta_r]],
                      [B[phi, delta_a], B[phi, delta_r]],
                      [B[psi, delta_a], B[psi, delta_r]]])    
    return A_lon, B_lon, A_lat, B_lat

def euler_state(x_quat):
    # convert state x with attitude represented by quaternion
    # to x_euler with attitude represented by Euler angles
    phi, theta, psi = quaternion_to_euler(x_quat[6:10])
    x_euler = np.zeros((12,1))
    x_euler[:6] = x_quat[:6]
    x_euler[6:9] = np.array([[phi, theta, psi]]).T
    x_euler[9:] = x_quat[10:]
    
    return x_euler

def quaternion_state(x_euler):
    # convert state x_euler with attitude represented by Euler angles
    # to x_quat with attitude represented by quaternions

    ##### TODO #####
    phi = x_euler.item(6)
    theta= x_euler.item(7)
    psi = x_euler.item(8)
    e = euler_to_quaternion(phi, theta, psi)
    
    x_quat = np.zeros((13,1))
    x_quat[:6] = x_euler[:6]
    x_quat[6:10] = e
    x_quat[10:] = x_euler[9:]
    return x_quat

def f_euler(mav, x_euler, delta):
    # return 12x1 dynamics (as if state were Euler state)
    # compute f at euler_state, f_euler will be f, except for the attitude states

    # need to correct attitude states by multiplying f by
    # partial of quaternion_to_euler(quat) with respect to quat
    # compute partial quaternion_to_euler(quat) with respect to quat
    # dEuler/dt = dEuler/dquat * dquat/dt
    x_quat = quaternion_state(x_euler)
    mav._state = x_quat
    mav._update_velocity_data()
    ##### TODO #####
    f_euler_ = np.zeros((12,1))

    return f_euler_

def df_dx(mav, x_euler, delta):
    # take partial of f_euler with respect to x_euler
    eps = 0.01  # deviation

    ##### TODO #####
    # Jacobian of f wrt x
    A = np.zeros((12, 12))
    f_at_x = f_euler(mav, x_euler, delta)
    for i in range(0, 12):
        x_eps = np.copy(x_euler)
        x_eps[i][0] += eps
        f_at_x_eps = f_euler(mav, x_eps, delta)
        df_dxi = (f_at_x_eps - f_at_x) / eps
        A[:, i] = df_dxi[:, 0]
    return A


def df_du(mav, x_euler, delta):
    # take partial of f_euler with respect to input
    eps = 0.01  # deviation

    ##### TODO #####
    B = np.zeros((12, 4))  # Jacobian of f wrt u
    f_at_u = f_euler(mav, x_euler, delta)
    for i in range(0, 4):
        delta_eps = np.copy(delta)
        delta_eps[i] += eps
        f_at_u_eps = f_euler(mav, x_euler, delta_eps)
        df_dui = (f_at_u_eps - f_at_u) / eps
        B[:, i] = df_dui[:, 0]
    return B


def dT_dVa(mav, Va, delta_t):
    # returns the derivative of motor thrust with respect to Va
    eps = 0.01

    ##### TODO #####
    dT_dVa = 0
    return dT_dVa

def dT_ddelta_t(mav, Va, delta_t):
    # returns the derivative of motor thrust with respect to delta_t
    eps = 0.01

    ##### TODO #####
    dT_ddelta_t = 0
    return dT_ddelta_t
