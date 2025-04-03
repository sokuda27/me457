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
end_time = 60

roll_hold = []
course_hold = []
sideslip_hold = []

while sim_time < end_time:

    # -------autopilot commands-------------
    commands.airspeed_command = Va_command.step(sim_time)
    commands.course_command = course_command.step(sim_time)
    # commands.altitude_command = altitude_command.square(sim_time)

    # -------autopilot-------------
    estimated_state = mav.true_state  # uses true states in the control
    delta, commanded_state = autopilot.update(commands, estimated_state)

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    mav.update(delta, current_wind)  # propagate the MAV dynamics

    euler = euler_state(mav._state[0:13])

    roll_hold.append(euler.item(6))
    course_hold.append(np.arctan(euler.item(3)/euler.item(4)))
    sideslip_hold.append(np.arctan(euler.item(4)/np.sqrt((euler.item(5)**2 + euler.item(3)**2))))

    # -------increment time-------------
    sim_time += SIM.ts_simulation
    # time.sleep(0.002) # slow down the simulation for visualization

    print("running")

fig, axs = plt.subplots(1, 3, constrained_layout=True, figsize=(15,5)) 
axs[0].plot(roll_hold)
axs[0].set_xlabel('Time (s)')
axs[0].set_ylabel('roll hold')

axs[1].plot(course_hold)
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('course hold')

# axs[2].plot(sideslip_hold)
# axs[2].set_xlabel('Time (s)')
# axs[2].set_ylabel('sideslip')

fig.suptitle('Position graphs', fontsize=16)
plt.show()