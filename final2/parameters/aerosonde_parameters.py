import numpy as np
from tools.rotations import euler_to_quaternion
# THIS VEHICLE IS CALLED PENGUIN B
######################################################################################
                #   Initial Conditions
######################################################################################
#   Initial conditions for MAV
north0 = 0.  # initial north position
east0 = 0.  # initial east position
down0 = -100.0  # initial down position
u0 = 25.  # initial velocity along body x-axis
v0 = 0.  # initial velocity along body y-axis
w0 = 0.  # initial velocity along body z-axis
phi0 = 0.  # initial roll angle
theta0 = 0.  # initial pitch angle
psi0 = 0.0  # initial yaw angle
p0 = 0  # initial roll rate
q0 = 0  # initial pitch rate
r0 = 0  # initial yaw rate
Va0 = np.sqrt(u0**2+v0**2+w0**2)
#   Quaternion State
e = euler_to_quaternion(phi0, theta0, psi0)
e0 = e.item(0)
e1 = e.item(1)
e2 = e.item(2)
e3 = e.item(3)


######################################################################################
                #   Physical Parameters
######################################################################################
mass = 24.0               # kg (typical with payload)
Jx = 1.2                  # kg*m^2
Jy = 1.7                  # kg*m^2
Jz = 2.7                  # kg*m^2
Jxz = 0.15                # kg*m^2 (small cross-product inertia)
S_wing = 1.86             # m^2 (wing area)
b = 3.3                   # m (wingspan)
c = S_wing / b            # m (mean aerodynamic chord, approx)
S_prop = 0.25             # m^2 (estimated)
rho = 1.225               # kg/m^3 (standard air density at sea level)
e = 0.85                  # Oswald efficiency factor
AR = b**2 / S_wing        # Aspect Ratio
gravity = 9.81            # m/s^2

######################################################################################
                #   Longitudinal Coefficients
######################################################################################
C_L_0 = 0.2
C_D_0 = 0.04
C_m_0 = 0.01
C_L_alpha = 5.5
C_D_alpha = 0.1
C_m_alpha = -2.8
C_L_q = 7.5
C_D_q = 0.0
C_m_q = -36.0
C_L_delta_e = 0.12
C_D_delta_e = 0.01
C_m_delta_e = -1.0
M = 50.0                  # Sigma parameter for blending (from McLain)
alpha0 = 0.47             # Angle of attack at zero lift for nonlinear blending
epsilon = 0.16            # Induced drag factor
C_D_p = 0.043             # Parasite drag coefficient


######################################################################################
                #   Lateral Coefficients
######################################################################################
C_Y_0 = 0.0
C_ell_0 = 0.0
C_n_0 = 0.0
C_Y_beta = -0.9
C_ell_beta = -0.15
C_n_beta = 0.07
C_Y_p = 0.0
C_ell_p = -0.5
C_n_p = 0.07
C_Y_r = 0.0
C_ell_r = 0.25
C_n_r = -0.09
C_Y_delta_a = 0.07
C_ell_delta_a = 0.17
C_n_delta_a = -0.01
C_Y_delta_r = 0.18
C_ell_delta_r = 0.002
C_n_delta_r = -0.07

######################################################################################
                #   Propeller thrust / torque parameters (see addendum by McLain)
######################################################################################

D_prop = 0.635           # Propeller diameter (25 inches converted to meters)
KV_rpm_per_volt = 150.0  # Motor speed constant (RPM per Volt)
KV = (1. / KV_rpm_per_volt) * 60. / (2. * np.pi)  # Back-emf constant in V-s/rad
KQ = KV                  # Motor torque constant, N-m/A
R_motor = 0.04           # Motor resistance (Ohms)
i0 = 1.5                 # No-load current (Amps)
ncells = 12.0            # 12S battery
V_max = 3.7 * ncells     # Max voltage (Volts)

# Coeffiecients from prop_data fit
C_Q2 = -0.01664
C_Q1 = 0.004970
C_Q0 = 0.005230
C_T2 = -0.1079
C_T1 = -0.06044
C_T0 = 0.09357

######################################################################################
                #   Calculation Variables
######################################################################################
#   gamma parameters pulled from page 36 (dynamics)
gamma = Jx * Jz - (Jxz**2)
gamma1 = (Jxz * (Jx - Jy + Jz)) / gamma
gamma2 = (Jz * (Jz - Jy) + (Jxz**2)) / gamma
gamma3 = Jz / gamma
gamma4 = Jxz / gamma
gamma5 = (Jz - Jx) / Jy
gamma6 = Jxz / Jy
gamma7 = ((Jx - Jy) * Jx + (Jxz**2)) / gamma
gamma8 = Jx / gamma

#   C values defines on pag 62
C_p_0         = gamma3 * C_ell_0      + gamma4 * C_n_0
C_p_beta      = gamma3 * C_ell_beta   + gamma4 * C_n_beta
C_p_p         = gamma3 * C_ell_p      + gamma4 * C_n_p
C_p_r         = gamma3 * C_ell_r      + gamma4 * C_n_r
C_p_delta_a    = gamma3 * C_ell_delta_a + gamma4 * C_n_delta_a
C_p_delta_r    = gamma3 * C_ell_delta_r + gamma4 * C_n_delta_r
C_r_0         = gamma4 * C_ell_0      + gamma8 * C_n_0
C_r_beta      = gamma4 * C_ell_beta   + gamma8 * C_n_beta
C_r_p         = gamma4 * C_ell_p      + gamma8 * C_n_p
C_r_r         = gamma4 * C_ell_r      + gamma8 * C_n_r
C_r_delta_a    = gamma4 * C_ell_delta_a + gamma8 * C_n_delta_a
C_r_delta_r    = gamma4 * C_ell_delta_r + gamma8 * C_n_delta_r