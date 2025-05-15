from plotter.plotter import Plotter
import random
import numpy as np
from PyQt5.QtWidgets import QApplication, QWidget  
import sys

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