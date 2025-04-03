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
from controllers.autopilot_lqr import Autopilot
from models.compute_models import euler_state
# from controllers.autopilot_lqr import Autopilot
import time

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)
autopilot = Autopilot(SIM.ts_simulation)

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
course_command = Signals(dc_offset=np.radians(180),
                         amplitude=np.radians(45),
                         start_time=5.0,
                         frequency=0.015)

sim_time = SIM.start_time
end_time = 30

roll_hold = []
pitch_hold = []
yaw_hold = []
course_hold = []
course_command_plot = []


while sim_time < end_time:

    # -------autopilot commands-------------
    # commands.airspeed_command = Va_command.step(sim_time)
    commands.course_command = course_command.step(sim_time)
    course_command_plot.append(course_command.step(sim_time))
    # commands.altitude_command = altitude_command.square(sim_time)

    # -------autopilot-------------
    estimated_state = mav.true_state  # uses true states in the control
    delta, commanded_state = autopilot.update(commands, estimated_state)

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    mav.update(delta, current_wind)  # propagate the MAV dynamics

    roll_hold.append(mav.true_state.phi)
    pitch_hold.append(mav.true_state.theta)
    yaw_hold.append(mav.true_state.psi)
    course_hold.append(mav.true_state.chi)

    # -------increment time-------------
    sim_time += SIM.ts_simulation
    # time.sleep(0.002) # slow down the simulation for visualization

fig, axs = plt.subplots(1, 3, constrained_layout=True, figsize=(15,5))

axs[0].plot(course_hold)
axs[0].plot(course_command_plot)
axs[0].set_xlabel('Time (s)')
axs[0].set_ylabel('roll hold')
axs[0].legend()

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