import numpy as np
import matplotlib.pyplot as plt

from models.mav_dynamics_control import MavDynamics
from models.wind_simulation import WindSimulation
from message_types.msg_delta import MsgDelta
import parameters.simulation_parameters as SIM
from models.trim import compute_trim
from tools.rotations import euler_to_quaternion
from tools.signals import Signals
from models.mav_dynamics_control import MavDynamics
from models.wind_simulation import WindSimulation
from controllers.autopilot import Autopilot
from models.compute_models import compute_model
from models.compute_models import compute_ss_model
# from controllers.autopilot_lqr import Autopilot
from scipy import signal
from numpy import linalg as LA
import models.model_coef as TF
import time

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)
autopilot = Autopilot(SIM.ts_simulation)

trim_state_in, trim_input_in = compute_trim(mav, 25, 0)
compute_model(mav, trim_state_in, trim_input_in)

# A_lon, B_lon, A_lat, B_lat = compute_ss_model(mav, 25, 0)

# sys_lat = signal.StateSpace(TF.A_lat, TF.B_lat, np.zeros((12,12)), np.zeros((12,4)))

[s_lat, v_lat] = LA.eig(TF.A_lat)

print(s_lat)

