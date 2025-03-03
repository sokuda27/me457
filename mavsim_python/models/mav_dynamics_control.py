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
        wind_body = quaternion_to_rotation(self._state[6:10]) @ steady_state 
        # add the gust 
        wind_body += gust
        # convert total wind to world frame
        self._wind = self._state[3:5] - wind_body

        # velocity vector relative to the airmass ([ur , vr, wr]= ?)
        [ur, vr, wr] = self._wind
        # compute airspeed (self._Va = ?)
        self._Va = np.sqrt(ur**2 + vr**2 + wr**2)
        # compute angle of attack (self._alpha = ?)
        self._alpha = np.arctan(wr/ur)
        # compute sideslip angle (self._beta = ?)
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
        p,q,r = self._state[10:12]
        
        # compute gravitational forces ([fg_x, fg_y, fg_z])
        fg_x = -MAV.mass*9.8*np.sin(theta)
        fg_y = MAV.mass*9.8*np.cos(theta)*np.sin(phi)
        fg_z = MAV.mass*9.8*np.cos(theta)*np.cos(phi)

        # compute Lift and Drag coefficients (CL, CD)
        CL = MAV.C_L_0 + MAV.C_L_alpha*self._alpha
        CD = MAV.C_D_0 + MAV.C_D_alpha*self._alpha
        
        # compute Lift and Drag Forces (F_lift, F_drag)
        F_lift = 0.5*MAV.rho*self._Va**2*MAV.S_wing*CL
        F_drag = 0.5*MAV.rho*self._Va**2*MAV.S_wing*CD
        
        # propeller thrust and torque
        thrust_prop, torque_prop = self._motor_thrust_torque(self._Va, delta.throttle)

        # compute longitudinal forces in body frame (fx, fz)
        fx = -np.cos(self._alpha)*(F_drag) + np.sin(self._alpha)*F_lift
        fz = -np.sin(self._alpha)*(F_drag) - np.cos(self._alpha)*F_lift
        
        
        # compute lateral forces in body frame (fy)
        fy = 
        
        # compute logitudinal torque in body frame (My)

        # compute lateral torques in body frame (Mx, Mz)

        forces_moments = np.array([[0, 0, 0, 0, 0, 0]]).T
        return forces_moments

    def _motor_thrust_torque(self, Va, delta_t):
        # compute thrust and torque due to propeller
        ##### TODO #####
        # map delta_t throttle command(0 to 1) into motor input voltage
        # v_in =

        # Angular speed of propeller (omega_p = ?)

        # thrust and torque due to propeller
        thrust_prop = 0
        torque_prop = 0

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
