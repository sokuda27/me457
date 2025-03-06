import numpy as np
from models.mav_dynamics import MavDynamics
from models.wind_simulation import WindSimulation
import parameters.simulation_parameters as SIM
from models.trim import compute_trim

wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)

Va = 25.
gamma = 0.*np.pi/180.
trim_state, trim_input = compute_trim(mav, Va, gamma)
mav._state = trim_state  # set the initial state of the mav to the trim state
delta = trim_input  # set input to constant constant trim input

sim_time = SIM.start_time
end_time = 60 

while sim_time < end_time:
    current_wind = wind.update() 

# intg.__doc__
# plt.figure()
# plt.title('pn')
# p_north = [arr[0] for arr in x_history]
# p_east = [arr[1] for arr in x_history]
# p_down = [arr[2] for arr in x_history]
# plt.plot(t_history, p_north)
# #plt.legend(["pn", "pe", "pd"])
# plt.show()

# plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot3D(p_north, p_east, p_down)
# ax.set_title('Path of aircraft')
# ax.set_xlabel('north', fontsize=10)
# ax.set_ylabel('east', fontsize=10)
# ax.set_zlabel('down', fontsize=10)
# plt.show()
