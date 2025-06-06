"""
        1/5/2023 - David L. Christiansen
        7/13/2023 - RWB
"""
import os, sys
# insert parent directory at beginning of python search path
from pathlib import Path
sys.path.insert(0,os.fspath(Path(__file__).parents[2]))
import numpy as np
from models.mav_dynamics_control import MavDynamics
from models.compute_models_old import compute_ss_model, compute_tf_model, euler_state, quaternion_state, f_euler, df_dx, df_du, dT_dVa, dT_ddelta_t
import parameters.simulation_parameters as SIM
from message_types.msg_delta import MsgDelta

mav = MavDynamics(SIM.ts_simulation)

trim_state = np.array([[0.000000, -0.000000, -100.000000, 24.968743, 0.000000, 1.249755, 0.999687, 0.000000, 0.025003, 0.000000, 0.000000, 0.000000, 0.000000]]).T
trim_input = MsgDelta(elevator=-0.124778,
                          aileron=0.001836,
                          rudder=-0.000303,
                          throttle=0.676752)
mav._state = trim_state

A_lon, B_lon, A_lat, B_lat = compute_ss_model(mav, trim_state, trim_input)

print("A_lon:\n" , A_lon , "\n")
print("B_lon:\n" , B_lon, "\n")
print("A_lat:\n" , A_lat, "\n")
print("B_lat:\n" , B_lat, "\n")

Va_trim, alpha_trim, theta_trim, a_phi1, a_phi2, a_theta1, a_theta2, a_theta3, \
    a_V1, a_V2, a_V3 = compute_tf_model(mav, trim_state, trim_input)

print("    Va_trim: " , Va_trim, "\n")
print(" alpha_trim: " , alpha_trim, "\n")
print(" theta_trim: " , theta_trim, "\n")
print("     a_phi1: " , a_phi1, "\n")
print("     a_phi2: " , a_phi2, "\n")
print("   a_theta1: " , a_theta1, "\n")
print("   a_theta2: " , a_theta2, "\n")
print("   a_theta3: " , a_theta3, "\n")
print("       a_V1: " , a_V1, "\n")
print("       a_V2: " , a_V2, "\n")
print("       a_V3: " , a_V3, "\n")

x_euler_ = euler_state(trim_state)
print("x_euler_: \n" , x_euler_, "\n")

x_quat_ = quaternion_state(x_euler_)
print("x_quat_: \n" , x_quat_, "\n")

f_euler_ = f_euler(mav, x_euler_, trim_input)
print("f_euler_: \n" , f_euler_, "\n")

A_ = df_dx(mav, x_euler_, trim_input)
print("A_: \n" , A_, "\n")

B_ = df_du(mav, x_euler_, trim_input)
print("B_: \n" , B_, "\n")

dT_dVa_ =  dT_dVa(mav, mav._Va, trim_input.throttle)
print("dT_dVa_: " , dT_dVa_, "\n")

dT_ddelta_t_ =  dT_ddelta_t(mav, mav._Va, trim_input.throttle)
print("dT_ddelta_t_: " , dT_ddelta_t_, "\n\n\n")


mav._state = np.array([[ 6.19506532e+01],
 [ 2.22940203e+01],
 [-1.10837551e+02],
 [ 2.73465947e+01],
 [ 6.19628233e-01],
 [ 1.42257772e+00],
 [ 9.38688796e-01],
 [ 2.47421558e-01],
 [ 6.56821468e-02],
 [ 2.30936730e-01],
 [ 4.98772167e-03],
 [ 1.68736005e-01],
 [ 1.71797313e-01]])

mav._Va = 22.4
new_state = mav._state
delta = MsgDelta()
delta.elevator = -0.2
delta.aileron = 0.0
delta.rudder = 0.005
delta.throttle = 0.5
print(" Second Case \n\n")

x_euler_ = euler_state(new_state)
print("x_euler_: \n" , x_euler_, "\n")

x_quat_ = quaternion_state(x_euler_)
print("x_quat_: \n" , x_quat_, "\n")

f_euler_ = f_euler(mav, x_euler_, delta)
print("f_euler_: \n" , f_euler_, "\n")

A_ = df_dx(mav, x_euler_, delta)
print("A_: \n" , A_, "\n")

B_ = df_du(mav, x_euler_, delta)
print("B_: \n" , B_, "\n")

dT_dVa_ =  dT_dVa(mav, mav._Va, delta.throttle)
print("dT_dVa_: " , dT_dVa_, "\n")

dT_ddelta_t_ =  dT_ddelta_t(mav, mav._Va, delta.throttle)
print("dT_ddelta_t_: " , dT_ddelta_t_, "\n")

# A_lon:
#  [[-0.20676658  0.50039026 -1.21983882 -9.79511927 -0.        ]
#  [-0.56064206 -4.46393561 24.37105023 -0.53938541 -0.        ]
#  [ 0.19993539 -3.99297865 -5.29473836  0.         -0.        ]
#  [ 0.          0.          0.99997406  0.         -0.        ]
#  [ 0.04999035 -0.9987497  -0.         24.99958361  0.        ]] 

# B_lon:
#  [[ -0.13840016   8.20722086]
#  [ -2.58618345   0.        ]
#  [-36.11239041   0.        ]
#  [  0.           0.        ]
#  [ -0.          -0.        ]] 

# A_lat:
#  [[-7.76772629e-01  1.24975500e+00 -2.49687430e+01  9.79757127e+00
#    0.00000000e+00]
#  [-3.86671935e+00 -2.26288510e+01  1.09050409e+01  0.00000000e+00
#    0.00000000e+00]
#  [ 7.83077145e-01 -1.15091678e-01 -1.22765475e+00  0.00000000e+00
#    0.00000000e+00]
#  [ 0.00000000e+00  9.99999666e-01  5.00528958e-02  0.00000000e+00
#    0.00000000e+00]
#  [ 0.00000000e+00 -1.67051761e-08  1.00125153e+00  0.00000000e+00
#    0.00000000e+00]] 

# B_lat:
#  [[  1.48617191   3.76496884]
#  [130.88368125  -1.79637441]
#  [  5.01173513 -24.88134191]
#  [  0.           0.        ]
#  [  0.           0.        ]] 

#     Va_trim:  25.000000291201477 

#  alpha_trim:  0.05001104395214544 

#  theta_trim:  0.05001119284259128 

#      a_phi1:  22.62885095683996 

#      a_phi2:  130.88368124853207 

#    a_theta1:  5.294738359662443 

#    a_theta2:  99.94742395724161 

#    a_theta3:  -36.11239040790846 

#        a_V1:  0.2817100013494505 

#        a_V2:  8.20722086381344 

#        a_V3:  9.809999999999892 

# x_euler_: 
#  [[ 0.00000000e+00]
#  [-0.00000000e+00]
#  [-1.00000000e+02]
#  [ 2.49687430e+01]
#  [ 0.00000000e+00]
#  [ 1.24975500e+00]
#  [ 0.00000000e+00]
#  [ 5.00111928e-02]
#  [ 0.00000000e+00]
#  [ 0.00000000e+00]
#  [ 0.00000000e+00]
#  [ 0.00000000e+00]] 

# x_quat_: 
#  [[ 0.00000000e+00]
#  [-0.00000000e+00]
#  [-1.00000000e+02]
#  [ 2.49687430e+01]
#  [ 0.00000000e+00]
#  [ 1.24975500e+00]
#  [ 9.99687376e-01]
#  [ 0.00000000e+00]
#  [ 2.50029906e-02]
#  [ 0.00000000e+00]
#  [ 0.00000000e+00]
#  [ 0.00000000e+00]
#  [ 0.00000000e+00]] 

# f_euler_: 
#  [[ 2.50000003e+01]
#  [ 0.00000000e+00]
#  [-3.72226119e-06]
#  [-5.12768098e-04]
#  [ 1.58782607e-03]
#  [ 9.98742522e-03]
#  [ 0.00000000e+00]
#  [ 0.00000000e+00]
#  [ 0.00000000e+00]
#  [-4.98388573e-05]
#  [-1.47473075e-06]
#  [ 2.51707628e-04]] 

# A_: 
#  [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00  9.98749701e-01
#    0.00000000e+00  4.99903481e-02 -3.12375834e-04 -1.25002682e-01
#   -1.24998960e-01  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
#    1.00000000e+00  0.00000000e+00 -1.24973417e+00  0.00000000e+00
#    2.49995836e+01  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00 -4.99903481e-02
#    0.00000000e+00  9.98749701e-01 -6.24091015e-03 -2.49995836e+01
#    1.38750773e-14  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00 -2.06766579e-01
#   -3.63078202e-05  5.00390258e-01  0.00000000e+00 -9.79511927e+00
#    0.00000000e+00  0.00000000e+00 -1.21983882e+00  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.26892670e-04
#   -7.76772629e-01  6.37546449e-06  9.79757127e+00  0.00000000e+00
#    0.00000000e+00  1.24975500e+00  0.00000000e+00 -2.49687430e+01]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00 -5.60642055e-01
#   -1.56603950e-04 -4.46393561e+00 -4.89882646e-02 -5.39385406e-01
#    1.29236899e-13  0.00000000e+00  2.43710502e+01  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
#    0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
#    0.00000000e+00  9.99999666e-01  0.00000000e+00  5.00528958e-02]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
#    0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
#    0.00000000e+00 -2.43928799e-05  9.99974056e-01 -2.56443516e-05]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
#    0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
#    0.00000000e+00 -1.67051761e-08  0.00000000e+00  1.00125153e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.43288661e-01
#   -3.86671935e+00  7.19895526e-03  0.00000000e+00  0.00000000e+00
#    0.00000000e+00 -2.26288510e+01  0.00000000e+00  1.09050409e+01]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.99935388e-01
#   -2.35956915e-11 -3.99297865e+00  0.00000000e+00  0.00000000e+00
#    0.00000000e+00 -1.06079295e-03 -5.29473836e+00  1.06079295e-03]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  9.82820776e-03
#    7.83077145e-01  4.93778316e-04  0.00000000e+00  0.00000000e+00
#    0.00000000e+00 -1.15091678e-01  0.00000000e+00 -1.22765475e+00]] 

# B_: 
#  [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [-1.38400158e-01  0.00000000e+00  0.00000000e+00  8.20722086e+00]
#  [ 0.00000000e+00  1.48617191e+00  3.76496884e+00  1.30104261e-16]
#  [-2.58618345e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [ 0.00000000e+00  1.30883681e+02 -1.79637441e+00 -5.31386382e+00]
#  [-3.61123904e+01  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [ 0.00000000e+00  5.01173513e+00 -2.48813419e+01 -3.63723254e-01]] 

# dT_dVa_:  -2.3521982958853327 

# dT_ddelta_t_:  90.27942950194785 



#  Second Case 


# x_euler_: 
#  [[ 6.19506532e+01]
#  [ 2.22940203e+01]
#  [-1.10837551e+02]
#  [ 2.73465947e+01]
#  [ 6.19628233e-01]
#  [ 1.42257772e+00]
#  [ 5.17674540e-01]
#  [ 9.03286236e-03]
#  [ 4.84851312e-01]
#  [ 4.98772167e-03]
#  [ 1.68736005e-01]
#  [ 1.71797313e-01]] 

# x_quat_: 
#  [[ 6.19506532e+01]
#  [ 2.22940203e+01]
#  [-1.10837551e+02]
#  [ 2.73465947e+01]
#  [ 6.19628233e-01]
#  [ 1.42257772e+00]
#  [ 9.38688796e-01]
#  [ 2.47421558e-01]
#  [ 6.56821468e-02]
#  [ 2.30936730e-01]
#  [ 4.98772167e-03]
#  [ 1.68736005e-01]
#  [ 1.71797313e-01]] 

# f_euler_: 
#  [[ 2.42832387e+01]
#  [ 1.26051301e+01]
#  [ 1.29573271e+00]
#  [-1.31778614e+00]
#  [-3.41507931e-01]
#  [ 1.24861387e+00]
#  [ 7.06527544e-03]
#  [ 6.16229478e-02]
#  [ 2.32773972e-01]
#  [ 2.20041055e-01]
#  [ 2.05034334e+00]
#  [ 2.12128940e-01]] 

# A_: 
#  [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00  8.84708165e-01
#   -4.01053084e-01  2.37587641e-01  7.17285818e-01  1.02534355e+00
#   -1.27263352e+01  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  4.66057800e-01
#    7.70901599e-01 -4.34166848e-01 -1.36496690e+00  5.40143494e-01
#    2.42198088e+01  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00 -9.03273952e-03
#    4.94840529e-01  8.68936857e-01 -1.73242216e-01 -2.73654375e+01
#    8.88178420e-14  0.00000000e+00  0.00000000e+00  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00 -2.06437917e-01
#    1.67758476e-01  3.92286478e-01  2.22044605e-14 -9.80899325e+00
#    2.22044605e-14  0.00000000e+00 -1.38851585e+00  6.19628233e-01]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00 -1.89364602e-01
#   -8.52344326e-01  4.07064881e-03  8.49985677e+00 -6.81211843e-02
#    0.00000000e+00  1.42257772e+00  0.00000000e+00 -2.73465947e+01]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00 -4.41741693e-01
#   -2.46890492e-02 -4.89162969e+00 -4.89692568e+00 -1.19620371e-01
#    1.33226763e-13 -6.19628233e-01  2.66918142e+01  0.00000000e+00]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
#    0.00000000e+00  0.00000000e+00  4.57332278e-04  2.32785596e-01
#    6.84717270e-06  9.99803542e-01  4.65945661e-03  7.52232695e-03]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
#    0.00000000e+00  0.00000000e+00 -2.33045276e-01 -1.45951852e-04
#    2.88111478e-05  9.72629246e-05  8.68868663e-01 -4.94693211e-01]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
#    0.00000000e+00  0.00000000e+00  6.04210285e-02  3.21737153e-03
#   -4.96668711e-05 -1.81663545e-04  4.95056161e-01  8.68704383e-01]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.22933671e-01
#   -4.23813655e+00  6.41593064e-03  0.00000000e+00  0.00000000e+00
#    0.00000000e+00 -2.47721978e+01 -1.32477696e-01  1.18171038e+01]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  4.12603256e-01
#    4.22266864e-03 -4.36582619e+00  0.00000000e+00  0.00000000e+00
#    0.00000000e+00  1.39345132e-01 -5.80103825e+00  4.16161389e-02]
#  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.02393302e-02
#    8.59089374e-01  5.34512684e-04  0.00000000e+00  0.00000000e+00
#    0.00000000e+00 -1.54489161e-01 -2.17077302e-02 -1.36554366e+00]] 

# B_: 
#  [[  0.           0.           0.           0.        ]
#  [  0.           0.           0.           0.        ]
#  [  0.           0.           0.           0.        ]
#  [ -0.16004176   0.           0.           5.51743963]
#  [  0.           1.78398623   4.51943179   0.        ]
#  [ -3.10474937   0.           0.           0.        ]
#  [  0.           0.           0.           0.        ]
#  [  0.           0.           0.           0.        ]
#  [  0.           0.           0.           0.        ]
#  [  0.         157.1114916   -2.15635028  -4.51605243]
#  [-43.34896045   0.           0.           0.        ]
#  [  0.           6.01603786 -29.86731962  -0.30911467]] 

# dT_dVa_:  -2.3737948095355677 

# dT_ddelta_t_:  60.69183591548715 
