import numpy as np
import matplotlib.pyplot as plt

from models.mav_dynamics_control import MavDynamics
from models.wind_simulation import WindSimulation
from message_types.msg_delta import MsgDelta
import parameters.simulation_parameters as SIM
from models.trim import compute_trim
from tools.rotations import euler_to_quaternion

wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)

# Va = 25.
# gamma = 0.*np.pi/180.
# trim_state, trim_input = compute_trim(mav, Va, gamma)
# mav._state = trim_state  # set the initial state of the mav to the trim state
# delta = trim_input  # set input to constant constant trim input

phi0 = 0; theta0 = 0; psi0 = 0

e = euler_to_quaternion(phi0, theta0, psi0)
e0 = e.item(0)
e1 = e.item(1)
e2 = e.item(2)
e3 = e.item(3)


mav._state = np.array([[0], [0], [-10], [25], [0], [0], [e0], [e1], [e2], [e3], [0], [0], [0], [0], [0]])
delta = MsgDelta()

sim_time = SIM.start_time
end_time = 60 

# with trim
# delta_e_trim = delta.elevator
# delta_a_trim = delta.aileron
# delta_r_trim = delta.rudder

# no input
delta.elevator = 0
delta.aileron = 0
delta.rudder = 0

pnorth_hist = [mav._state[0]]
peast_hist = [mav._state[1]]
pdown_hist = [mav._state[2]]
time = [sim_time]

while sim_time < end_time:
    current_wind = wind.update()
    mav.update(delta, current_wind)

    pnorth_hist.append(mav._state.item(0))
    peast_hist.append(mav._state.item(1))
    pdown_hist.append(mav._state.item(2))

    sim_time += SIM.ts_simulation
    time.append(sim_time)

fig, axs = plt.subplots(1, 3, constrained_layout=True, figsize=(15,5)) 
axs[0].plot(time[2:], pnorth_hist[2:])
axs[0].set_xlabel('Time (s)')
axs[0].set_ylabel('pn')

axs[1].plot(time[2:], peast_hist[2:])
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('pe')

axs[2].plot(time[2:], pdown_hist[2:])
axs[2].set_xlabel('Time (s)')
axs[2].set_ylabel('pd')

fig.suptitle('Position graphs', fontsize=16)
plt.show()

plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(pnorth_hist[2:], peast_hist[2:], pdown_hist[2:])
ax.set_title('Path of aircraft')
ax.set_xlabel('north', fontsize=10)
ax.set_ylabel('east', fontsize=10)
ax.set_zlabel('down', fontsize=10)
plt.show()