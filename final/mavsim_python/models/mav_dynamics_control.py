"""
mavDynamics 
    - this file implements the dynamic equations of motion for MAV
    - use unit quaternion for the attitude state
    
mavsim_python
    - Beard & McLain, PUP, 2012
    - Update history:  
        2/24/2020 - RWB
"""
import numpy as np
import math
from models.mav_dynamics import MavDynamics as MavDynamicsForces
# load message types
from message_types.msg_state import MsgState
from message_types.msg_delta import MsgDelta
import parameters.aerosonde_parameters as MAV
from tools.rotations import quaternion_to_rotation, quaternion_to_euler


class MavDynamics(MavDynamicsForces):
    def __init__(self, Ts):
        super().__init__(Ts)
        # store wind data for fast recall since it is used at various points in simulation
        self._wind = np.array([[0.], [0.], [0.]])  # wind in NED frame in meters/sec
        # store forces to avoid recalculation in the sensors function
        self._forces = np.array([[0.], [0.], [0.]])
        self._Va = MAV.u0
        self._alpha = 0
        self._beta = 0
        # update velocity data and forces and moments
        self._update_velocity_data()
        self._forces_moments(delta=MsgDelta())
        # update the message class for the true state
        self._update_true_state()


    ###################################
    # public functions
    def update(self, delta, wind):
        '''
            Integrate the differential equations defining dynamics, update sensors
            delta = (delta_a, delta_e, delta_r, delta_t) are the control inputs
            wind is the wind vector in inertial coordinates
            Ts is the time step between function calls.
        '''
        # get forces and moments acting on rigid bod
        forces_moments = self._forces_moments(delta)
        super()._rk4_step(forces_moments)
        # update the airspeed, angle of attack, and side slip angles using new state
        self._update_velocity_data(wind)
        # update the message class for the true state
        self._update_true_state()

    ###################################
    # private functions
    def _update_velocity_data(self, wind=np.zeros((6,1))):
        steady_state = wind[0:3]
        gust = wind[3:6]

        ##### TODO #####
        # convert steady-state wind vector from world to body frame
        wind_body = quaternion_to_rotation(self._state[6:10]) @ steady_state + gust
        # add the gust - done
        
        # convert total wind to world frame
        
        self._wind = self._state[3:6] - wind_body

        # velocity vector relative to the airmass ([ur , vr, wr]= ?)
        [ur, vr, wr] = self._wind
        # compute airspeed (self._Va = ?)
        self._Va = np.sqrt(ur**2 + vr**2 + wr**2)
        # compute angle of attack (self._alpha = ?)
        if ur == 0:
            self._alpha = np.pi/2
        else:
            self._alpha = np.arctan(wr/ur)
        # compute sideslip angle (self._beta = ?)
        if self._Va == 0:
            self._beta = 0
        else:
            self._beta = np.arcsin(vr/self._Va)
        
    def _forces_moments(self, delta):
        """
        return the forces on the UAV based on the state, wind, and control surfaces
        :param delta: np.matrix(delta_a, delta_e, delta_r, delta_t)
        :return: Forces and Moments on the UAV np.matrix(Fx, Fy, Fz, Ml, Mn, Mm)
        """
        ##### TODO ######
        # extract states (phi, theta, psi, p, q, r)
        
        phi, theta, psi = quaternion_to_euler(self._state[6:10])
        p,q,r = self._state[10:13]
        
        # compute gravitational forces ([fg_x, fg_y, fg_z])
        fg = quaternion_to_rotation(self._state[6:10]).T @ np.array([[0], [0], [MAV.mass * MAV.gravity]])
        fg_x = fg[0]
        fg_y = fg[1]
        fg_z = fg[2]

        # compute Lift and Drag coefficients (CL, CD)
        CL = MAV.C_L_0 + MAV.C_L_alpha*self._alpha
        CD = MAV.C_D_0 + MAV.C_D_alpha*self._alpha
        
        # compute Lift and Drag Forces (F_lift, F_drag)
        c = MAV.c
        C_L_q = MAV.C_L_q
        C_L_delta_e = MAV.C_L_delta_e
        C_D_delta_e = MAV.C_D_delta_e
        C_D_q = MAV.C_D_q
        delta_e = delta.elevator
        F_lift = 0.5*MAV.rho*self._Va**2*MAV.S_wing*(CL+C_L_q*c*q/(2*self._Va)+C_L_delta_e*delta_e)
        F_drag = 0.5*MAV.rho*self._Va**2*MAV.S_wing*(CD+C_D_q*c*q/(2*self._Va)+C_D_delta_e*delta_e)
        
        # propeller thrust and torque
        thrust_prop, torque_prop = self._motor_thrust_torque(self._Va, delta.throttle)

        # compute longitudinal forces in body frame (fx, fz)
        fx = -np.cos(self._alpha)*(F_drag) + np.sin(self._alpha)*F_lift + fg_x + thrust_prop
        fz = -np.sin(self._alpha)*(F_drag) - np.cos(self._alpha)*F_lift + fg_z
        
        
        # compute lateral forces in body frame (fy)
        b = MAV.b
        C_Y_0 = MAV.C_Y_0
        C_Y_beta = MAV.C_Y_beta
        C_Y_p = MAV.C_Y_p
        C_Y_r = MAV.C_Y_r
        C_Y_delta_a = MAV.C_Y_delta_a
        C_Y_delta_r = MAV.C_Y_delta_r
        delta_a = delta.aileron
        delta_r = delta.rudder
        fy = 0.5*MAV.rho*self._Va**2*MAV.S_wing*(C_Y_0 + C_Y_beta*self._beta + (b/(2*self._Va))*(C_Y_p*p + C_Y_r*r) + C_Y_delta_a*delta_a + C_Y_delta_r*delta_r) + fg_y
        
        # compute logitudinal torque in body frame (My) m
        C_m_0 = MAV.C_m_0
        C_m_alpha = MAV.C_m_alpha
        C_m_q = MAV.C_m_q
        C_m_delta_e = MAV.C_m_delta_e
        My = 0.5*MAV.rho*self._Va**2*MAV.S_wing*c*(C_m_0 + C_m_alpha*self._alpha + C_m_q*c*q/(2*self._Va) + C_m_delta_e*delta_e)

        # compute lateral torques in body frame (Mx, Mz) (l,n)
        C_ell_0 = MAV.C_ell_0
        C_ell_beta = MAV.C_ell_beta
        C_ell_p = MAV.C_ell_p
        C_ell_r = MAV.C_ell_r
        C_ell_delta_a = MAV.C_ell_delta_a
        C_ell_delta_r = MAV.C_ell_delta_r
        C_n_0 = MAV.C_n_0
        C_n_beta = MAV.C_n_beta
        C_n_p = MAV.C_n_p
        C_n_r = MAV.C_n_r
        C_n_delta_a = MAV.C_n_delta_a
        C_n_delta_r = MAV.C_n_delta_r
        Mx = 0.5*MAV.rho*self._Va**2*MAV.S_wing*b*(C_ell_0 + C_ell_beta*self._beta + (b/(2*self._Va))*(C_ell_p*p + C_ell_r*r) + C_ell_delta_a*delta_a + C_ell_delta_r*delta_r) - torque_prop
        Mz = 0.5*MAV.rho*self._Va**2*MAV.S_wing*b*(C_n_0 + C_n_beta*self._beta + (b/(2*self._Va))*(C_n_p*p + C_n_r*r) + C_n_delta_a*delta_a + C_n_delta_r*delta_r)

        forces_moments = np.array([[fx, fy, fz, Mx, My, Mz]]).T #saira note: check order!!!
        return forces_moments

    def _motor_thrust_torque(self, Va, delta_t):
        # compute thrust and torque due to propeller
        ##### TODO #####
        # map delta_t throttle command(0 to 1) into motor input voltage
        rho = MAV.rho
        v_max = MAV.V_max
        d_prop = MAV.D_prop
        kq = MAV.KQ
        R = MAV.R_motor
        io = MAV.i0
        cq2 = MAV.C_Q2
        cq1 = MAV.C_Q1
        cqo = MAV.C_Q0
        ct2 = MAV.C_T2
        ct1 = MAV.C_T1
        ct0 = MAV.C_T0
        V_in = delta_t * v_max

        # Angular speed of propeller
        a = rho * (d_prop ** 5) * cqo / ((2. * np.pi) ** 2)
        b = ((rho * (d_prop ** 4) * cq1 * self._Va) / (2 * np.pi)) + (kq **2 / R)
        c = (rho * (d_prop ** 3) * cq2 * self._Va ** 2) + (kq * io) - (kq * V_in / R)
        Omega_p = (-b + (b ** 2 - 4 * a * c)**0.5) / (2.0 * a)
        J = 2 * np.pi * Va / (Omega_p * d_prop)
        # thrust and torque due to propeller
        ct_func = ct2 * J ** 2 + ct1 * J + ct0
        cq_func = cq2 * J ** 2 + cq1 * J + cqo
        n = Omega_p / (2.0 * np.pi)
        # thrust and torque due to propeller
        #thrust_prop = 0.5 * MAV.rho * MAV.S_prop * d_prop * ((R * delta_t)**2 - Va**2)
        thrust_prop = MAV.rho * n**2 * np.power(d_prop, 4) * ct_func
        #torque_prop = kq * ((1 / R) * (V_in - kv * Omega_p) - io)
        torque_prop = MAV.rho * n**2 * np.power(d_prop, 5) * cq_func
        #torque_prop = MAV.rho * MAV.C_Q0 * MAV.D_prop**5 / (4 * math.pi**2) * sigma_p**2 + MAV.rho * MAV.C_Q1 * Va * MAV.D_prop**4 / (2 * math.pi) * sigma_p + MAV.rho * MAV.C_Q2 * MAV.D_prop**3 * Va**2
        return thrust_prop, torque_prop

    def _update_true_state(self):
        # rewrite this function because we now have more information
        phi, theta, psi = quaternion_to_euler(self._state[6:10])
        pdot = quaternion_to_rotation(self._state[6:10]) @ self._state[3:6]
        self.true_state.north = self._state.item(0)
        self.true_state.east = self._state.item(1)
        self.true_state.altitude = -self._state.item(2)
        self.true_state.Va = self._Va
        self.true_state.alpha = self._alpha
        self.true_state.beta = self._beta
        self.true_state.phi = phi
        self.true_state.theta = theta
        self.true_state.psi = psi
        self.true_state.Vg = np.linalg.norm(pdot)
        self.true_state.gamma = np.arcsin(pdot.item(2) / self.true_state.Vg)
        self.true_state.chi = np.arctan2(pdot.item(1), pdot.item(0))
        self.true_state.p = self._state.item(10)
        self.true_state.q = self._state.item(11)
        self.true_state.r = self._state.item(12)
        self.true_state.wn = self._wind.item(0)
        self.true_state.we = self._wind.item(1)
        self.true_state.bx = 0
        self.true_state.by = 0
        self.true_state.bz = 0
        self.true_state.camera_az = 0
        self.true_state.camera_el = 0
