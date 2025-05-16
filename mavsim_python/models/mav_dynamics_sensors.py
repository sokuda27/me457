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
from message_types.msg_sensors import MsgSensors
import parameters.aerosonde_parameters as MAV
import parameters.sensor_parameters as SENSOR
from models.mav_dynamics_control import MavDynamics as MavDynamicsNoSensors
from tools.rotations import quaternion_to_rotation, quaternion_to_euler, euler_to_rotation
from math import cos
from math import sin

class MavDynamics(MavDynamicsNoSensors):
    def __init__(self, Ts):
        super().__init__(Ts)
        # initialize the sensors message
        self._sensors = MsgSensors()
        # random walk parameters for GPS
        self._gps_eta_n = 0.
        self._gps_eta_e = 0.
        self._gps_eta_h = 0.
        # timer so that gps only updates every ts_gps seconds
        self._t_gps = 999.  # large value ensures gps updates at initial time.

    def sensors(self):
        "Return value of sensors on MAV: gyros, accels, absolute_pressure, dynamic_pressure, GPS"
        state = self.true_state
        phi, theta, psi = quaternion_to_euler(self._state[6:10])
        p, q, r = self._state[10], self._state[11], self._state[12]
        u, v, w = self._state[3], self._state[4], self._state[5]
        fx, fy, fz = self._forces.item(0), self._forces.item(1), self._forces.item(2)
        h = -self._state[2]
        m = MAV.mass
        g = MAV.gravity
        B_abs_pres = 0
        B_diff_pres = 0

        # simulate rate gyros (rad/s)
        self._sensors.gyro_x = float(p + SENSOR.gyro_x_bias + np.random.normal(0, SENSOR.gyro_sigma))
        self._sensors.gyro_y = float(q + SENSOR.gyro_y_bias + np.random.normal(0, SENSOR.gyro_sigma))
        self._sensors.gyro_z = float(r + SENSOR.gyro_z_bias + np.random.normal(0, SENSOR.gyro_sigma))

        # simulate accelerometers (m/s^2)
        self._sensors.accel_x = float((fx/m) + g * np.sin(theta) + np.random.normal(0, SENSOR.accel_sigma))
        self._sensors.accel_y = float((fy/m) - g * np.cos(theta) * np.sin(phi) + np.random.normal(0, SENSOR.accel_sigma))
        self._sensors.accel_z = float((fz/m) - g * np.cos(theta) * np.cos(phi) + np.random.normal(0, SENSOR.accel_sigma))

        # simulate magnetometers
        a = np.radians(12.5)
        b = np.radians(66)
        m_i = np.array([
            np.cos(a) * np.cos(b),
            np.cos(a) * np.sin(b),
            np.sin(a)])

        m_b = quaternion_to_rotation(self._state[6:10]) @ m_i
        self._sensors.mag_x = float(m_b.item(0) + np.random.normal(0, SENSOR.mag_sigma))
        self._sensors.mag_y = float(m_b.item(1) + np.random.normal(0, SENSOR.mag_sigma))
        self._sensors.mag_z = float(m_b.item(2) + np.random.normal(0, SENSOR.mag_sigma))

        # simulate pressure sensors
        self._sensors.abs_pressure = float(MAV.rho * MAV.gravity * h + B_abs_pres + np.random.normal(0, SENSOR.abs_pres_sigma))
        self._sensors.diff_pressure = float(MAV.rho * (self._Va**2)/2 + B_diff_pres + np.random.normal(0, SENSOR.diff_pres_sigma))

        # simulate GPS sensor
        if self._t_gps >= SENSOR.ts_gps:
            self._gps_eta_n = np.exp(-SENSOR.gps_k * SENSOR.ts_gps) * self._gps_eta_n + np.random.normal(0, SENSOR.gps_n_sigma)
            self._gps_eta_e = np.exp(-SENSOR.gps_k * SENSOR.ts_gps) * self._gps_eta_e + np.random.normal(0, SENSOR.gps_e_sigma)
            self._gps_eta_h = np.exp(-SENSOR.gps_k * SENSOR.ts_gps) * self._gps_eta_h + np.random.normal(0, SENSOR.gps_h_sigma)
            self._sensors.gps_n = float(state.north + self._gps_eta_n)
            self._sensors.gps_e = float(state.east + self._gps_eta_e)
            self._sensors.gps_h = float(-self._state.item(2) + self._gps_eta_h)
            self._sensors.gps_Vg = float(np.sqrt((self._Va * np.sin(psi) + state.wn)**2 + (self._Va * np.cos(psi) + state.we)**2) + np.random.normal(0, SENSOR.gps_Vg_sigma))
            self._sensors.gps_course = float(np.arctan2(self._Va * np.sin(psi) + state.wn, self._Va * np.cos(psi) + state.we) + np.random.normal(0, SENSOR.gps_course_sigma))
            self._t_gps = 0.0
        else:
            self._t_gps += self._ts_simulation
        return self._sensors

    def external_set_state(self, new_state):
        self._state = new_state

    def _update_true_state(self):
        # update the class structure for the true state:
        #   [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
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
        self.true_state.bx = SENSOR.gyro_x_bias
        self.true_state.by = SENSOR.gyro_y_bias
        self.true_state.bz = SENSOR.gyro_z_bias
        # self.true_state.camera_az = self._state.item(13)
        # self.true_state.camera_el = self._state.item(14)
