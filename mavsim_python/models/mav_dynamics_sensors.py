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
        pn = state.north
        pe = state.east
        h = -self._state.item(2)
        #pd = h
        p = self._state.item(10)
        q = self._state.item(11)
        r = self._state.item(12)
        phi = state.phi
        theta = state.theta
        psi = state.psi
        g = MAV.gravity
        fx = self._forces[0]
        fy = self._forces[1]
        fz = self._forces[2]
        Va = state.Va
        # simulate rate gyros(units are rad / sec)
        self._sensors.gyro_x = p + SENSOR.gyro_x_bias + np.random.normal(0, SENSOR.gyro_sigma)
        self._sensors.gyro_y = q + SENSOR.gyro_y_bias + np.random.normal(0, SENSOR.gyro_sigma)
        self._sensors.gyro_z = r + SENSOR.gyro_z_bias + np.random.normal(0, SENSOR.gyro_sigma)

        # simulate accelerometers(units of g)
        self._sensors.accel_x = (fx / MAV.mass + g * sin(theta) + np.random.normal(0, SENSOR.accel_sigma)).item(0)
        self._sensors.accel_y = (
                    fy / MAV.mass - g * cos(theta) * sin(phi) + np.random.normal(0, SENSOR.accel_sigma)).item(0)
        self._sensors.accel_z = (
                    fz / MAV.mass - g * cos(theta) * cos(phi) + np.random.normal(0, SENSOR.accel_sigma)).item(0)

        # simulate magnetometers
        # magnetic field in provo has magnetic declination of 12.5 degrees
        # and magnetic inclination of 66 degrees
        
        # simulate magnetometers
        # magnetic field in Provo: declination = 12.5 deg east, inclination = 66 deg downward
        inclination = np.radians(66.0)
        declination = np.radians(12.5)
        mag_strength = 1.0  # Assume normalized strength of 1 (Tesla units don't matter here)

        # magnetic field in inertial (north-east-down) frame
        mag_n = mag_strength * np.cos(inclination) * np.cos(declination)
        mag_e = mag_strength * np.cos(inclination) * np.sin(declination)
        mag_d = mag_strength * np.sin(inclination)

        mag_inertial = np.array([[mag_n], [mag_e], [mag_d]])

        # rotation from body to inertial
        R = quaternion_to_rotation(self._state[6:10])  # body-to-inertial rotation matrix
        R_transpose = R.T  # inertial-to-body frame

        mag_body = R_transpose @ mag_inertial  # rotate magnetic field into body frame

        # add optional small measurement noise
        self._sensors.mag_x = mag_body.item(0) + np.random.normal(0, SENSOR.mag_sigma)
        self._sensors.mag_y = mag_body.item(1) + np.random.normal(0, SENSOR.mag_sigma)
        self._sensors.mag_z = mag_body.item(2) + np.random.normal(0, SENSOR.mag_sigma)

        # simulate pressure sensors
        h_AGL = h
        self._sensors.abs_pressure = MAV.rho * MAV.gravity * h_AGL + np.random.normal(0,
                                                                                                         SENSOR.abs_pres_sigma)
        self._sensors.diff_pressure = MAV.rho * self._Va ** 2 / 2 + np.random.normal(0,
                                                                                                      SENSOR.diff_pres_sigma)
        
        # simulate GPS sensor
        if self._t_gps >= SENSOR.ts_gps:
            self._gps_eta_n = np.exp(-SENSOR.gps_k * SENSOR.ts_gps) * self._gps_eta_n + np.random.normal(0,
                                                                                                         SENSOR.gps_n_sigma)
            self._gps_eta_e = np.exp(-SENSOR.gps_k * SENSOR.ts_gps) * self._gps_eta_e + np.random.normal(0,
                                                                                                         SENSOR.gps_e_sigma)
            self._gps_eta_h = np.exp(-SENSOR.gps_k * SENSOR.ts_gps) * self._gps_eta_h + np.random.normal(0,
                                                                                                         SENSOR.gps_h_sigma)
            self._sensors.gps_n = pn + self._gps_eta_n
            self._sensors.gps_e = pe + self._gps_eta_e
            self._sensors.gps_h = -self._state.item(2) + self._gps_eta_h
            wn = state.wn
            we = state.we
            V1 = Va * np.sin(psi) + we
            V2 = Va * np.cos(psi) + wn
            self._sensors.gps_Vg = (np.sqrt(V1 ** 2 + V2 ** 2) + np.random.normal(0, SENSOR.gps_Vg_sigma)).item(0)
            self._sensors.gps_course = (np.arctan2(V1, V2) + np.random.normal(0, SENSOR.gps_course_sigma)).item(0)
            self._t_gps = 0.
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
        self.true_state.camera_az = self._state.item(13)
        self.true_state.camera_el = self._state.item(14)    
