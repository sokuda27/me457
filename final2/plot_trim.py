from plotter.plotter import Plotter
import random
import numpy as np
from PyQt5.QtWidgets import QApplication, QWidget  
import sys

from models.mav_dynamics_control import MavDynamics
from models.wind_simulation import WindSimulation
from message_types.msg_delta import MsgDelta
import parameters.simulation_parameters as SIM
from models.trim import compute_trim
from tools.rotations import euler_to_quaternion
from tools.rotations import quaternion_to_euler
from tools.signals import Signals
from models.mav_dynamics_sensors import MavDynamics
from models.wind_simulation import WindSimulation
from controllers.autopilot import Autopilot
from models.compute_models import compute_model
# from controllers.autopilot_lqr import Autopilot
from estimators.observer import Observer
import time

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

trim_state_in, trim_input_in = compute_trim(mav, 25, 0)
hello = compute_model(mav, trim_state_in, trim_input_in)

p = Plotter(QApplication(sys.argv), 1)
p.show_window()
p.create_plot_widget(plot_id="trim_speed", xlabel="time", ylabel="V_a",background_color='w')
p.create_data_set(plot_id="trim_speed",data_label="v_a",data_color=(0,255,0),data_thickness=20)
# p.set_plot_data("trim_speed","tree", #time# , #whatever va over time is#)
# p.add_data_points("fig","tree", np.linspace(50,60,50).tolist(), np.linspace(0,200,50).tolist())
for i in range(200):
    p.add_data_point("foo","truth",i,random.randint(0,20))
    p.add_data_point("fi","truth",i,random.randint(0,20))
    p.add_data_point("fi","estimate",i,random.randint(0,20))
    # time.sleep(0.1)
    p.update_plots()
    if i == 20:
        p.save_image()
p.hold_window_until_exit()