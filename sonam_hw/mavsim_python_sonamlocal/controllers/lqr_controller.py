import sys
from numpy import array, sin, cos, radians, concatenate, zeros, diag
from scipy.linalg import solve_continuous_are, inv
sys.path.append('..')
import parameters.control_parameters as AP
from tools.wrap import wrap
import models.mav_dynamics_control as MAV
from message_types.msg_state import MsgState
from message_types.msg_delta import MsgDelta

class lqr:
    def __init__(self, ts_control):
        self.Ts = ts_control
        # initialize integrators and delay variables
        self.integratorCourse = 0
        self.integratorAltitude = 0
        self.integratorAirspeed = 0
        self.errorCourseD1 = 0
        self.errorAltitudeD1 = 0
        self.errorAirspeedD1 = 0
        #compute LQR gains
        CrLat = array([[0, 0, 0, 0, 1.0]])
        AAlat = concatenate(())