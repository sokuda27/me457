# Parameters file for the Aerosonde

# Physical Parameters
mass = 11.0       # Mass in kg
g = 9.81          # Gravitational Acceleration in m/s^2
J_x = 0.824       # Moment of Inertia about x in kg*m^2
J_y = 1.135       # Moment of Inertia about y in kg*m^2
J_z = 1.759       # Moment of Inertia about z in kg*m^2
J_xy = 0          # Product of Inertia in kg*m^2 (Assumed to be 0 for symmetric body)
J_yz = 0          # Product of Inertia in kg*m^2 (Assumed to be 0 for symmetric body)
J_xz = 0          # Product of Inertia in kg*m^2
S = 0.55          # Wing Area in m^2
b = 2.9           # Wing Span in m
c = 0.19          # Mean Aerodynamic Chord in m
rho = 1.268       # Air Density in kg/m^3
e = 0.9           # Oswald Efficiency Factor

# Motor Parameters
V_max = 44.4      # Maximum Velocity in m/s
D_prop = 0.508    # Diameter of Propeller in m
K_V = 0.0659      # V*s/rad
K_Q = 0.0659      # in N-m
Omega = V_max/K_V # Propellor Speed in rad/s
R_motor = 0.01    # Resistance of motor in ohms
i_0 = 1.5         # No load current of motor in A
C_T0 = 0.09357    # Thrust Coefficient at hover
C_T1 = -0.06044   # Thrust Coefficient at forward flight
C_T2 = -0.1079    # Thrust Coefficient at high speeds
C_Q0 = 0.005230   # Torque Coefficient at hover
C_Q1 = 0.004970   # Torque Coefficient at forward flight
C_Q2 = -0.01664   # Torque Coefficient at high speeds

# Aerodynamic Longitudinal Coefficients
C_L_0 = 0.23
C_D_0 = 0.043
C_m_0 = 0.0135
C_L_alpha = 5.61
C_D_alpha = 0.030
C_m_alpha = -2.74
C_L_q = 7.95
C_D_q = 0
C_m_q = -38.21
C_L_delta_e = 0.13
C_D_delta_e = 0.0135
C_m_delta_e = -0.99
M = 50
alpha_0 = 0.47
C_D_p = 0.043

# Aerodynamic Lateral Coefficients
C_Y_0 = 0
C_l_0 = 0
C_n_0 = 0
C_Y_beta = -0.83
C_l_beta = -0.12
C_n_beta = 0.073
C_Y_p = 0
C_l_p = -0.51
C_n_p = -0.069
C_Y_r = 0
C_l_r = 0.25
C_n_r = -0.095
C_Y_delta_a = 0.075
C_l_delta_a = 0.17
C_n_delta_a = -0.011
C_Y_delta_r = 0.19
C_l_delta_r = 0.0024
C_n_delta_r = -0.069