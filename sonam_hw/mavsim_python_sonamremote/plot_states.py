import numpy as np
import matplotlib.pyplot as plt

from models.mav_dynamics_control import MavDynamics
from models.wind_simulation import WindSimulation
import parameters.simulation_parameters as SIM
from models.trim import compute_trim

wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)

Va = 25.
gamma = 0.*np.pi/180.
trim_state, trim_input = compute_trim(mav, Va, gamma)
mav._state = np.array([0,0,0,0,0,0,0,0,0,0,0,0]) # trim_state  # set the initial state of the mav to the trim state
delta = trim_input  # set input to constant constant trim input

print(trim_input)

sim_time = SIM.start_time
end_time = 60 

# with trim
delta_e_trim = delta.elevator
delta_a_trim = delta.aileron
delta_r_trim = delta.rudder

# no input
# delta.elevator = 0
# delta.aileron = 0
# delta.rudder = 0

pnorth_hist = [mav._state[0]]
peast_hist = [mav._state[1]]
pdown_hist = [mav._state[2]]
time = [0]

while sim_time < end_time:
    current_wind = wind.update()
    mav.update(delta, current_wind)

    time.append(sim_time)

    pnorth_hist.append(mav._state.item(0))
    peast_hist.append(mav._state.item(1))
    pdown_hist.append(mav._state.item(2))

    sim_time += SIM.ts_simulation

# plt.figure()
# plt.plot(time, pnorth_hist)
# plt.show()

plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(pnorth_hist, peast_hist, pdown_hist)
ax.set_title('Path of aircraft')
ax.set_xlabel('north', fontsize=10)
ax.set_ylabel('east', fontsize=10)
ax.set_zlabel('down', fontsize=10)
plt.show()