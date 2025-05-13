"""
observer
    - Beard & McLain, PUP, 2012
    - Last Update:
        3/2/2019 - RWB
        3/4/2024 - RWB
"""
import numpy as np
import parameters.control_parameters as CTRL
import parameters.sensor_parameters as SENSOR
import parameters.aerosonde_parameters as MAV
from tools.wrap import wrap
from message_types.msg_state import MsgState
from message_types.msg_sensors import MsgSensors
from estimators.filters import AlphaFilter, ExtendedKalmanFilterContinuousDiscrete

class Observer:
    def __init__(self, ts: float, initial_measurements: MsgSensors=MsgSensors()):
        self.Ts = ts  # sample rate of observer
        # initialized estimated state message
        self.estimated_state = MsgState()

        ##### TODO #####
        self.lpf_gyro_x = AlphaFilter(alpha=0.5, y0=initial_measurements.gyro_x)
        self.lpf_gyro_y = AlphaFilter(alpha=0.5, y0=initial_measurements.gyro_y)
        self.lpf_gyro_z = AlphaFilter(alpha=0.5, y0=initial_measurements.gyro_z)
        self.lpf_accel_x = AlphaFilter(alpha=0.5, y0=initial_measurements.accel_x)
        self.lpf_accel_y = AlphaFilter(alpha=0.5, y0=initial_measurements.accel_y)
        self.lpf_accel_z = AlphaFilter(alpha=0.5, y0=initial_measurements.accel_z)
        # use alpha filters to low pass filter absolute and differential pressure
        self.lpf_abs = AlphaFilter(alpha=0.9, y0=initial_measurements.abs_pressure)
        self.lpf_diff = AlphaFilter(alpha=0.5, y0=initial_measurements.diff_pressure)
        # ekf for phi and theta
        self.attitude_ekf = ExtendedKalmanFilterContinuousDiscrete(
            f=self.f_attitude,
            Q=np.diag([np.radians(0.5)**2, np.radians(0.5)**2]),  # (2x2)
            P0=np.diag([np.radians(5.0)**2, np.radians(5.0)**2]),  # (2x2)
            xhat0=np.array([[0.0], [0.0]]),  # (2x1)
            Qu=np.diag([
                SENSOR.gyro_sigma**2,
                SENSOR.gyro_sigma**2,
                SENSOR.gyro_sigma**2,
                SENSOR.abs_pres_sigma
            ]),  # (4x4)
            Ts=ts,
            N=5
            )
        # ekf for pn, pe, Vg, chi, wn, we, psi
        self.position_ekf = ExtendedKalmanFilterContinuousDiscrete(
            f=self.f_smooth, 
            Q=np.diag([(1.0)**2]*7), 
            P0=np.diag([(5.0)**2]*7), 
            xhat0=np.zeros((7, 1)), 
            Qu=np.diag([
                SENSOR.gyro_sigma**2, 
                SENSOR.gyro_sigma**2, 
                SENSOR.abs_pres_sigma,
                1.0, # guess for noise on roll
                1.0, # guess for noise on pitch
                ]),
            Ts=ts,
            N=10
            )
        self.R_accel = np.diag([
                SENSOR.accel_sigma**2, 
                SENSOR.accel_sigma**2, 
                SENSOR.accel_sigma**2
                ])
        self.R_pseudo = np.diag([
                0.01,  # pseudo measurement #1 ##### TODO #####
                0.01,  # pseudo measurement #2 ##### TODO #####
                ])
        self.R_gps = np.diag([
                    SENSOR.gps_n_sigma**2,  # y_gps_n
                    SENSOR.gps_e_sigma**2,  # y_gps_e
                    SENSOR.gps_Vg_sigma**2,  # y_gps_Vg
                    SENSOR.gps_course_sigma**2,  # y_gps_course
                    ])
        self.gps_n_old = 9999
        self.gps_e_old = 9999
        self.gps_Vg_old = 9999
        self.gps_course_old = 9999

    def update(self, measurement: MsgSensors) -> MsgState:
        ##### TODO #####
        # estimates for p, q, r are low pass filter of gyro minus bias estimate
        self.estimated_state.p = self.lpf_gyro_y.update(measurement.gyro_x)
        self.estimated_state.q = self.lpf_gyro_y.update(measurement.gyro_y)
        self.estimated_state.r = self.lpf_gyro_z.update(measurement.gyro_z)
        # invert sensor model to get     altitude and airspeed
        abs_pressure = self.lpf_abs.update(measurement.abs_pressure)
        diff_pressure = self.lpf_diff.update(measurement.diff_pressure)
        self.estimated_state.altitude = (abs_pressure) / (MAV.rho * MAV.gravity)
        self.estimated_state.Va = np.sqrt(2 * diff_pressure / MAV.rho)
        # estimate phi and theta with ekf
        u_attitude=np.array([
                [self.estimated_state.p],
                [self.estimated_state.q],
                [self.estimated_state.r],
                [self.estimated_state.Va],
                ])
        xhat_attitude, _ = self.attitude_ekf.propagate_model(u_attitude)
        y_accel=np.array([
                [measurement.accel_x],
                [measurement.accel_y],
                [measurement.accel_z],
                ])
        xhat_attitude, _=self.attitude_ekf.measurement_update(
            y=y_accel, 
            u=u_attitude,
            h=self.h_accel,
            R=self.R_accel)
        self.estimated_state.phi = xhat_attitude.item(0)
        self.estimated_state.theta = xhat_attitude.item(1)
        # estimate pn, pe, Vg, chi, wn, we, psi with ekf
        u_smooth = np.array([
                [self.estimated_state.q],
                [self.estimated_state.r],
                [self.estimated_state.Va],
                [self.estimated_state.phi],
                [self.estimated_state.theta],
                ])
        xhat_position, P_position=self.position_ekf.propagate_model(u_smooth)
        y_pseudo = np.array([[0.], [0.]])
        xhat_position, P_position=self.position_ekf.measurement_update(
            y=y_pseudo,
            u=u_smooth,
            h=self.h_pseudo,
            R=self.R_pseudo)
        # only update GPS when one of the signals changes
        if (measurement.gps_n != self.gps_n_old) \
            or (measurement.gps_e != self.gps_e_old) \
            or (measurement.gps_Vg != self.gps_Vg_old) \
            or (measurement.gps_course != self.gps_course_old):
            y_gps = np.array([
                    [measurement.gps_n],
                    [measurement.gps_e],
                    [measurement.gps_Vg],
                    [wrap(measurement.gps_course, xhat_position.item(3))],
                    ])
            xhat_position, P_position=self.position_ekf.measurement_update(
                y=y_gps,
                u=u_smooth,
                h=self.h_gps,
                R=self.R_gps)
            # update stored GPS signals
            self.gps_n_old = measurement.gps_n
            self.gps_e_old = measurement.gps_e
            self.gps_Vg_old = measurement.gps_Vg
            self.gps_course_old = measurement.gps_course
        self.estimated_state.north = xhat_position.item(0)
        self.estimated_state.east = xhat_position.item(1)
        self.estimated_state.Vg = xhat_position.item(2)
        self.estimated_state.chi = xhat_position.item(3)
        self.estimated_state.wn = xhat_position.item(4)
        self.estimated_state.we = xhat_position.item(5)
        self.estimated_state.psi = xhat_position.item(6)
        # not estimating these
        self.estimated_state.alpha = self.estimated_state.theta
        self.estimated_state.beta = 0.0
        self.estimated_state.bx = 0.0
        self.estimated_state.by = 0.0
        self.estimated_state.bz = 0.0
        return self.estimated_state

    def f_attitude(self, x: np.ndarray, u: np.ndarray) -> np.ndarray:
        '''
            system dynamics for propagation model: xdot = f(x, u)
                x = [phi, theta].T
                u = [p, q, r, Va].T
        '''
        ##### TODO #####
        phi = x.item(0)
        theta = x.item(1)
        p = u.item(0)
        q = u.item(1)
        r = u.item(2)
        phi_dot = p + q*np.sin(phi)*np.tan(theta) + r*np.cos(phi)*np.tan(theta)
        theta_dot = q*np.cos(phi) - r*np.sin(phi)
        xdot = np.array([[phi_dot, theta_dot]])
        return xdot

    def h_accel(self, x: np.ndarray, u: np.ndarray)->np.ndarray:
        '''
            measurement model y=h(x,u) for accelerometers
                x = [phi, theta].T
                u = [p, q, r, Va].T
        '''
        ##### TODO #####
        phi = x.item(0)
        theta = x.item(1)
        p = u.item(0)
        q = u.item(1)
        r = u.item(2)
        Va = u.item(3)
        accel_x = q*Va*np.sin(theta) + SENSOR.gravity*np.sin(theta)
        accel_y = r*Va*np.cos(theta) - p*Va*np.sin(theta) - SENSOR.gravity*np.cos(theta)*np.sin(phi)
        accel_z = -q*Va*np.cos(theta) - SENSOR.gravity*np.cos(theta)*np.cos(phi)
        y = np.array([[accel_x], [accel_y], [accel_z]])
        return y

    def f_smooth(self, x: np.ndarray, u: np.ndarray) -> np.ndarray:
        '''
            system dynamics for propagation model: xdot = f(x, u)
                x = [pn, pe, Vg, chi, wn, we, psi].T
                u = [p, q, r, Va, phi, theta].T
        '''
        ##### TODO #####   
        pn = x.item(0)
        pe = x.item(1)
        Vg = x.item(2)
        chi = x.item(3)
        wn = x.item(4)
        we = x.item(5)
        psi = x.item(6)

        q = u.item(0)
        r = u.item(1)
        Va = u.item(2)
        phi = u.item(3)
        theta = u.item(4)

        pndot = Vg * np.cos(chi)
        pedot = Vg * np.sin(chi)
        chidot = SENSOR.gravity/Vg*np.tan(phi)*np.cos(chi-psi)
        psidot = q*np.sin(phi)/np.cos(theta) + r*np.cos(phi)/np.cos(theta)
        Vgdot = ((Va*np.cos(psi) + wn) * (-Va*psidot*np.sin(psi)) + (Va*np.sin(psi) + we)*(Va*psidot*np.cos(psi)))/Vg

        xdot = np.array([[pndot], [pedot], [Vgdot], [chidot], [0.0], [0.0], [psidot]])
        return xdot

    def h_pseudo(self, x: np.ndarray, u: np.ndarray)->np.ndarray:
        '''
            measurement model measurement model for wind triangale pseudo measurement: y=y(x, u)
                x = [pn, pe, Vg, chi, wn, we, psi].T
                u = [q, r, Va, phi, theta].T
            returns
                y = [pn, pe, Vg, chi]
        '''
        ##### TODO #####       
        Vg = x.item(2)
        chi = x.item(3)
        wn = x.item(4)
        we = x.item(5)
        psi = x.item(6)
        Va = u.item(2)
        y = np.array([
            [Va*np.cos(psi) + wn - Vg*np.cos(chi)],
            [Va*np.sin(psi) + we - Vg*np.sin(chi)]
        ])
        return y

    def h_gps(self, x: np.ndarray, u: np.ndarray)->np.ndarray:
        '''
            measurement model for gps measurements: y=y(x, u)
                x = [pn, pe, Vg, chi, wn, we, psi].T
                u = [p, q, r, Va, phi, theta].T
            returns
                y = [pn, pe, Vg, chi]
        '''
        ##### TODO ##### 
        pn = x.item(0)
        pe = x.item(1)
        Vg = x.item(2)
        chi = x.item(3)        
        y = np.array([
            [pn],
            [pe],
            [Vg],
            [chi]
        ])
        return y




