from models.mav_dynamics_control import MavDynamics
from models.wind_simulation import WindSimulation
from message_types.msg_delta import MsgDelta
from models.compute_models import euler_state
import parameters.simulation_parameters as SIM
from models.trim import compute_trim
from tools.rotations import euler_to_quaternion
import numpy as np

wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)

Va = 25.
gamma = 0.*np.pi/180.
trim_state, trim_input = compute_trim(mav, Va, gamma)
mav._state = trim_state  # set the initial state of the mav to the trim state
delta = trim_input  # set input to constant constant trim input
x_euler = euler_state(mav._state)

# with trim
delta_e_trim = delta.elevator
delta_a_trim = delta.aileron
delta_r_trim = delta.rudder

x_lat = np.array([[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                 [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]]) @ x_euler

print(x_lat)    