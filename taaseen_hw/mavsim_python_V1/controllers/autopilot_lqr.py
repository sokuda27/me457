"""
autopilot block for mavsim_python
    - Beard & McLain, PUP, 2012
    - Last Update:
        2/10/22 - RWB
"""
import numpy as np
from numpy import array, sin, cos, radians, concatenate, zeros, diag
from scipy.linalg import solve_continuous_are, inv
import parameters.control_parameters as AP
from tools.wrap import wrap
import models.model_coef as M
from message_types.msg_state import MsgState
from message_types.msg_delta import MsgDelta

def saturate(input, low_limit, up_limit):
    if input <= low_limit:
        output = low_limit
    elif input >= up_limit:
        output = up_limit
    else:
        output = input
    return output


class Autopilot:
    def __init__(self, ts_control):
        self.Ts = ts_control
        # initialize integrators and delay variables
        self.integratorCourse = 0
        self.integratorAltitude = 0
        self.integratorAirspeed = 0
        self.errorCourseD1 = 0
        self.errorAltitudeD1 = 0
        self.errorAirspeedD1 = 0
        # compute LQR gains
        
        #### TODO ######
        CrLat = array([[0, 0, 0, 0, 0]])
        AAlat = concatenate((
                    concatenate((M.A_lat, zeros((5,1))), axis=1),
                    concatenate((CrLat, zeros((1,1))), axis=1)),
                    axis=0)
        BBlat = concatenate((M.B_lat, zeros((1,2))), axis=0)
        Qlat = diag([0, 0, 0, 0, 0, 0]) # v, p, r, phi, chi, intChi
        Rlat = diag([0, 0]) # a, r
        # Plat = solve_continuous_are(AAlat, BBlat, Qlat, Rlat)
        Plat = Plon = np.zeros((6,6))
        # self.Klat = inv(Rlat) @ BBlat.T @ Plat
        self.Klat = np.zeros((2,6))
        CrLon = array([[0, 0, 0, 0, 0], [0, 0, 0, 0, 0]])
        AAlon = concatenate((
                    concatenate((M.A_lon, zeros((5,2))), axis=1),
                    concatenate((CrLon, zeros((2,2))), axis=1)),
                    axis=0)
        BBlon = concatenate((M.B_lon, zeros((2, 2))), axis=0)
        Qlon = diag([0, 0, 0, 0, 0, 0, 0]) # u, w, q, theta, h, intH, intVa
        Rlon = diag([0, 0])  # e, t
        # Plon = solve_continuous_are(AAlon, BBlon, Qlon, Rlon)
        Plon = np.zeros((7,7))
        # self.Klon = inv(Rlon) @ BBlon.T @ Plon
        self.Klat = np.zeros((2,7))
        self.commanded_state = MsgState()

    def update(self, cmd, state):
        # lateral autopilot
        errorAirspeed = state.Va - cmd.airspeed_command
        chi_c = wrap(cmd.course_command, state.chi)
        errorCourse = saturate(state.chi - chi_c, -radians(15), radians(15))
        self.integatorCourse = self.integratorCourse + (self.Ts/2) * (errorCourse + self.errorCourseD1)
        self.errorCourseD1 = errorCourse
        xLat = array([[errorAirspeed * sin(state.beta)], # v
                    [state.p],
                    [state.r],
                    [state.phi],
                    [errorCourse],
                    [self.integratorCourse]])

        tmp = -self.Klat @ xLat
        delta_a = saturate(tmp.item(0), -radians(30), radians(30))
        delta_r = saturate(tmp.item(1), -radians(30), radians(30))
        
        # longitudinal autopilot
        altitude_c = saturate(cmd.altitude_command,
                    state.altitude - 0.2*AP.altitude_zone,
                    state.altitude + 0.2*AP.altitude_zone)

        errorAltitude = state.altitude - altitude_c
        self.integatorAltitude = self.integratorAltitude + (self.Ts /2) * (errorAltitude + self.errorAltitudeD1)
        
        self.errorAltitudeD1 = errorAltitude
        self.integatorAirspeed = self.integratorAirspeed + (self.Ts/2) * (errorAirspeed + self.errorAirspeedD1)

        self.errorAirspeedD1 = errorAirspeed
        xLon = array ([[errorAirspeed * cos(state.alpha)], # u
                    [errorAirspeed * sin(state.alpha)], # w
                    [state.q],
                    [state.theta],
                    [errorAltitude],
                    [self.integratorAltitude],
                    [self.integratorAirspeed]])

        tmp = -self.Klon @ xLon
        delta_e = saturate(tmp.item(0), -radians(30), radians(30))
        delta_t = saturate(tmp.item(1), 0.0, 1.0)

        # construct control outputs and commanded states
        delta = MsgDelta(elevator=0,
                         aileron=0,
                         rudder=0,
                         throttle=0)
        self.commanded_state.altitude = 0
        self.commanded_state.Va = 0
        self.commanded_state.phi = 0
        self.commanded_state.theta = 0
        self.commanded_state.chi = 0
        return delta, self.commanded_state