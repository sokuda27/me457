"""
autopilot block for mavsim_python
    - Beard & McLain, PUP, 2012
    - Last Update:
        2/10/22 - RWB
"""
import os, sys
# insert parent directory at beginning of python search path
from pathlib import Path
sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
import numpy as np
from numpy import array, sin, cos, radians, concatenate, zeros, diag
from scipy.linalg import solve_continuous_are, inv
import parameters.control_parameters as AP
from tools.wrap import wrap
import models.model_coef as M
from message_types.msg_state import MsgState
from message_types.msg_delta import MsgDelta
from numpy.linalg import matrix_rank
import parameters.simulation_parameters as SIM
from parameters.simulation_parameters import ts_simulation as Ts

sim_time = SIM.start_time

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
        CrLat = array([[0, 0, 0, 0, 1]])
        AAlat = concatenate((
                    concatenate((M.A_lat, zeros((5,1))), axis=1),
                    concatenate((CrLat, zeros((1,1))), axis=1)),
                    axis=0)
        BBlat = concatenate((M.B_lat, zeros((1,2))), axis=0)
        Qlat = diag([1, 1, 1, 10, 10, 100]) # v, p, r, phi, chi, intChi
        Rlat = diag([1, 1]) # a, r
        # Plat = solve_continuous_are(AAlat, BBlat, Qlat, Rlat)
        # Plat = Plon = np.zeros((6,6))
        # self.Klat = inv(Rlat) @ BBlat.T @ Plat
        # self.Klat = np.zeros((2,6))
        CrLon = array([[0, 0, 0, 1, 0], [1/AP.Va0, 1/AP.Va0, 0, 0, 0]])
        AAlon = concatenate((
                    concatenate((M.A_lon, zeros((5,2))), axis=1),
                    concatenate((CrLon, zeros((2,2))), axis=1)),
                    axis=0)
        BBlon = concatenate((M.B_lon, zeros((2, 2))), axis=0)
        Qlon = diag([10, 10, 0.1, 1, 10, 100, 1000]) # u, w, q, theta, h, intH, intVa
        Rlon = diag([1, 1])  # e, t
        # Plon = solve_continuous_are(AAlon, BBlon, Qlon, Rlon)
        # Plon = np.zeros((7,7))
        # self.Klon = inv(Rlon) @ BBlon.T @ Plon
        # self.Klon = np.zeros((2,7))
        self.commanded_state = MsgState()         

        n = AAlat.shape[0]
        ctrb_matrix = BBlat
        for i in range(1, n):
            ctrb_matrix = np.hstack((ctrb_matrix, np.linalg.matrix_power(AAlat, i) @ BBlat))
        print(matrix_rank(ctrb_matrix) == n)
    
        n = AAlon.shape[0]
        ctrb_matrix = BBlon
        for i in range(1, n):
            ctrb_matrix = np.hstack((ctrb_matrix, np.linalg.matrix_power(AAlon, i) @ BBlon))
        print(matrix_rank(ctrb_matrix) == n)

autopilot = Autopilot(SIM.ts_simulation)