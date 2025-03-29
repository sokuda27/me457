import numpy as np
from scipy.optimize import minimize
from tools.rotations import euler_to_quaternion, quaternion_to_euler
import parameters.aerosonde_parameters as MAV
from parameters.simulation_parameters import ts_simulation as Ts
from message_types.msg_delta import MsgDelta

def compute_model(mav, trim_state, trim_input):
    A_lon, B_lon, A_lat, B_lat = compute_ss_model(mav, trim_state, trim_input)
    Va_trim, alpha_trim, theta_trim, a_phi1, a_phi2, a_theta1, a_theta2, a_theta3, \
    a_V1, a_V2, a_V3 = compute_tf_model(mav, trim_state, trim_input)

    # Write transfer function gains to file
    with open('models/model_coef.py', 'w') as file:
        file.write('import numpy as np\n')
        file.write(f'Va_trim = {Va_trim}\n')
        file.write(f'alpha_trim = {alpha_trim}\n')
        file.write(f'theta_trim = {theta_trim}\n')
        file.write(f'a_phi1 = {a_phi1}\n')
        file.write(f'a_phi2 = {a_phi2}\n')
        file.write(f'a_theta1 = {a_theta1}\n')
        file.write(f'a_theta2 = {a_theta2}\n')
        file.write(f'a_theta3 = {a_theta3}\n')
        file.write(f'a_V1 = {a_V1}\n')
        file.write(f'a_V2 = {a_V2}\n')
        file.write(f'a_V3 = {a_V3}\n')
        file.write(f'Ts = {Ts}\n')

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
    B_lon = np.array([[B[u,delta_e], B[u,delta_t]],
                      [B[w, delta_e], B[w,delta_t]],
                      [B[q,delta_e], B[q,delta_t]],
                      [B[theta,delta_e], B[theta,delta_t]],
                      [B[h,delta_e], B[h,delta_t]]])

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

def f_euler(mav, x_euler, delta):
    # Convert Euler state to quaternion state
    x_quat = quaternion_state(x_euler)
    mav._state = np.copy(x_quat)  # Ensure the state is correctly set in the MAV dynamics
    
    # Update velocity and forces/moments based on the new state
    mav._update_velocity_data()
    forces_moments = mav._forces_moments(delta)

    # Compute the state derivative using the MAV dynamics model
    x_dot_quat = mav._f(mav._state, forces_moments)

    # Convert quaternion state derivative back to Euler state derivative
    f_euler_ = euler_state(x_dot_quat)  

    return f_euler_

def df_dx(mav, x_euler, delta):
    eps = 0.001
    A = np.zeros((12, 12))
    f_at_x = f_euler(mav, x_euler, delta)
    for i in range(0,12):
        x_eps = np.copy(x_euler)
        x_eps[i] += eps
        f_at_x_eps = f_euler(mav, x_eps, delta)
        df_dxi = (f_at_x_eps - f_at_x) / eps
        A[:, i] = df_dxi [:]
    return A

def df_du(mav, x_euler, delta):
    eps = 0.01
    B = np.zeros((12, 4))  # 12 states, 4 control inputs (aileron, elevator, rudder, throttle)
    f_at_x = f_euler(mav, x_euler, delta)

    for i in range(4):
        delta_eps = MsgDelta(
            aileron=delta.aileron,
            elevator=delta.elevator,
            rudder=delta.rudder,
            throttle=delta.throttle,
        )

        # Apply perturbation to the i-th control input
        if i == 0:
            delta_eps.aileron += eps
        elif i == 1:
            delta_eps.elevator += eps
        elif i == 2:
            delta_eps.rudder += eps
        elif i == 3:
            delta_eps.throttle += eps

        f_at_x_eps = f_euler(mav, x_euler, delta_eps)
        B[:, i] = ((f_at_x_eps - f_at_x) / eps)

    return B

def dT_dVa(mav, Va, delta_t):
    """ Computes the derivative of thrust with respect to airspeed (Va) numerically """
    eps = 0.01  # Small perturbation

    if hasattr(mav, "_motor_thrust_torque"):
        T1, _ = mav._motor_thrust_torque(Va, delta_t)  # Extract only thrust
        T2, _ = mav._motor_thrust_torque(Va + eps, delta_t)
    else:
        raise AttributeError("MavDynamics object has no method '_motor_thrust_torque'. Check the correct function name.")

    return (T2 - T1) / eps  # Numerical derivative

def dT_ddelta_t(mav, Va, delta_t):
    """ Computes the derivative of thrust with respect to throttle (delta_t) numerically """
    eps = 0.01  # Small perturbation

    if hasattr(mav, "_motor_thrust_torque"):
        T1, _ = mav._motor_thrust_torque(Va, delta_t)  # Extract only thrust
        T2, _ = mav._motor_thrust_torque(Va, delta_t + eps)
    else:
        raise AttributeError("MavDynamics object has no method '_motor_thrust_torque'. Check the correct function name.")

    return (T2 - T1) / eps  # Numerical derivative

def euler_state(x_quat):
    euler = np.zeros(12)
    euler[:6] = x_quat[:6].flatten()  # Ensure 1D shape
    euler[6:9] = quaternion_to_euler(x_quat[6:10])
    euler[9:] = x_quat[10:].flatten()  # Ensure 1D shape
    return euler

def quaternion_state(x_euler):
    phi = x_euler.item(6)
    theta = x_euler.item(7)
    psi = x_euler.item(8)
    e = euler_to_quaternion(phi, theta, psi).reshape((4,1))  # Ensure correct shape

    x_quat = np.zeros((13,1))  # Ensure it's a column vector
    x_quat[:6] = x_euler[:6].reshape((6,1))  # Fix shape issue
    x_quat[6:10] = e
    x_quat[10:] = x_euler[9:].reshape((3,1))  # Fix shape issue
    return x_quat
