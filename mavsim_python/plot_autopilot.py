import numpy as np
import matplotlib.pyplot as plt

from models.mav_dynamics_control import MavDynamics
from models.wind_simulation import WindSimulation
from message_types.msg_delta import MsgDelta
import parameters.simulation_parameters as SIM
from models.trim import compute_trim
from tools.rotations import euler_to_quaternion
from tools.rotations import quaternion_to_euler
from tools.signals import Signals
from models.mav_dynamics_control import MavDynamics
from models.wind_simulation import WindSimulation
from controllers.autopilot import Autopilot
from models.compute_models import compute_model
from estimators.observer import Observer
import time

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)

# Initialize MAV to nominal state
mav._state = np.zeros((13, 1))
mav._Va = 25.0  # desired airspeed
mav._alpha = 0.0
mav._beta = 0.0
mav._wind = np.zeros((6, 1))
mav._forces = np.zeros((3, 1))

autopilot = Autopilot(SIM.ts_simulation)
observer = Observer(SIM.ts_simulation)

# autopilot commands
from message_types.msg_autopilot import MsgAutopilot
commands = MsgAutopilot()
Va_command = Signals(dc_offset=25.0,
                     amplitude=3.0,
                     start_time=2.0,
                     frequency=0.01)
altitude_command = Signals(dc_offset=100.0,
                           amplitude=20.0,
                           start_time=0.0,
                           frequency=0.02)
course_command = Signals(dc_offset=np.radians(45),
                         amplitude=np.radians(45),
                         start_time=5.0,
                         frequency=0.015)

sim_time = SIM.start_time
end_time = 30

roll_hold = []
roll_command_plot = []
# pitch_hold = []
# yaw_hold = []
course_hold = []
course_command_plot = []
alt_hold = []
alt_command_plot = []

# trim_state = np.array([[0.000000, -0.000000, -100.000000, 24.968743, 0.000000, 1.249755, 0.999687, 0.000000, 0.025003, 0.000000, 0.000000, 0.000000, 0.000000]]).T
# trim_input = MsgDelta(elevator=-0.124778,
#                           aileron=0.001836,
#                           rudder=-0.000303,
#                           throttle=0.676752)

trim_state_in, trim_input_in = compute_trim(mav, 25, 0)
hello = compute_model(mav, trim_state_in, trim_input_in)

while sim_time < end_time:

    # -------autopilot commands-------------
    commands.airspeed_command = Va_command.step(sim_time)
    commands.course_command = course_command.step(sim_time)
    commands.altitude_command = altitude_command.step(sim_time)
    course_command_plot.append(course_command.step(sim_time))
    alt_command_plot.append(commands.altitude_command)

    # -------autopilot-------------
    measurements = mav._state   
    estimated_state = mav.true_state  # uses true states in the control
    # print(quaternion_to_euler(measurements))
    delta, commanded_state = autopilot.update(commands, estimated_state)
    roll_command_plot.append(commanded_state.phi)

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    mav.update(delta, current_wind)  # propagate the MAV dynamics

    roll_hold.append(mav.true_state.phi)
    # pitch_hold.append(mav.true_state.theta)
    # yaw_hold.append(mav.true_state.psi)
    course_hold.append(mav.true_state.chi)
    alt_hold.append(mav.true_state.altitude)

    # -------increment time-------------
    sim_time += SIM.ts_simulation
    # time.sleep(0.002) # slow down the simulation for visualization

fig, axs = plt.subplots(1, 3, constrained_layout=True, figsize=(15,5))

axs[0].plot(course_hold)
axs[0].plot(course_command_plot)
axs[0].set_xlabel('Time (s)')
axs[0].set_ylabel('course hold')
axs[0].legend()

axs[1].plot(roll_hold)
axs[1].plot(roll_command_plot)
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('roll hold')
axs[1].legend()

axs[2].plot(alt_hold)
axs[2].plot(alt_command_plot)
axs[2].set_xlabel('Time (s)')
axs[2].set_ylabel('alt hold')
axs[2].legend()

# axs[0].plot(roll_hold)
# axs[0].set_xlabel('Time (s)')
# axs[0].set_ylabel('roll hold')

# axs[1].plot(pitch_hold)
# axs[1].set_xlabel('Time (s)')
# axs[1].set_ylabel('pitch hold')

# axs[2].plot(yaw_hold)
# axs[2].set_xlabel('Time (s)')
# axs[2].set_ylabel('yaw hold')

fig.suptitle('Position graphs', fontsize=16)
plt.show()