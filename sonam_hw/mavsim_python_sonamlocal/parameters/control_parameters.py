import numpy as np
import os, sys
# insert parent directory at beginning of python search path
from pathlib import Path
sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
import models.model_coef as TF
import parameters.aerosonde_parameters as MAV


#### TODO #####
gravity = MAV.gravity  # gravity constant
Va0 = TF.Va_trim
rho = MAV.rho # density of air
sigma = 0.005  # low pass filter gain for derivative

#----------roll loop-------------
# get transfer function data for delta_a to phi
wn_roll = 20 # error in roll when aileron saturates
zeta_roll = 0.707 # damping ratio for roll attitude loop
roll_kp = wn_roll**2/TF.a_phi2
roll_kd = (2.0 * zeta_roll * wn_roll - TF.a_phi1) / TF.a_phi2
# roll_kp = 0
# roll_kd = 0

#----------course loop-------------
wn_course = wn_roll / 20.0 # bandwidth between roll and course loops
zeta_course = 1.0 # damping ratio for course hold loop
course_kp = 2 * zeta_course * wn_course * Va0 / gravity
course_ki = wn_course**2 * Va0 / gravity

#----------yaw damper-------------
yaw_damper_p_wo = 0.4 # washout filter cutoff frequency
yaw_damper_kr = 0.2 # gain for yaw damper

#----------pitch loop-------------
wn_pitch = 24.0 # nat freq for pitch
zeta_pitch = 0.707 # damping ratio for pitch attitude loop
pitch_kp = (wn_pitch**2 - TF.a_theta2) / TF.a_theta3
pitch_kd = (2.0 * zeta_pitch * wn_pitch - TF.a_theta1) / TF.a_theta3
K_theta_DC = pitch_kp * TF.a_theta3 / (TF.a_theta2 + pitch_kp * TF.a_theta3)
# pitch_kp = 0
# pitch_kd = 0
# K_theta_DC = 1

#----------altitude loop-------------
wn_altitude =wn_pitch / 30.0 # bandwitdth separation between pitch and altitude
zeta_altitude = 1.0 # damping ratio for altitude hold loop
altitude_kp = 2.0 * zeta_altitude * wn_altitude / K_theta_DC / Va0
altitude_ki = wn_altitude**2 / K_theta_DC / Va0
altitude_zone = 10.0

#---------airspeed hold using throttle---------------
wn_airspeed_throttle = 3.0 # nat freq airspeed hold loop
zeta_airspeed_throttle = 2 # damping ratio for throttle loop
airspeed_throttle_kp = (2.0 * zeta_airspeed_throttle * wn_airspeed_throttle - TF.a_V1) / TF.a_V2
airspeed_throttle_ki = wn_airspeed_throttle**2 / TF.a_V2
# airspeed_throttle_kp = 0
# airspeed_throttle_ki = 0
