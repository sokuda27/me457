from models.mav_dynamics import MavDynamics
from models.wind_simulation import WindSimulation

mav = MavDynamics()
wind_sim = WindSimulation()
rk4 = mav.update()

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
