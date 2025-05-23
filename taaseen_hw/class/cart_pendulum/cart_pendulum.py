import numpy as np
import control as ct
import matplotlib.pyplot as plt

# Parameters
M = 0.5
m = 0.2
b = 0.1
I = 0.006
g = 9.8
l = 0.3

# State space
den = I*(M+m)+M*m*l**2
A = np.array([
    [0,      1,              0,           0],
    [0, -(I+m*l**2)*b/den,  (m**2*g*l**2)/den,  0],
    [0,      0,              0,           1],
    [0, -(m*l*b)/den,       m*g*l*(M+m)/den,  0]
    ])
B = np.array([
    [0],
    [(I+m*l**2)/den],
    [0],
    [m*l/den]
    ])
C = np.array([
    [1, 0, 0, 0],
    [0, 0, 1, 0]
    ])
D = np.array([
    [0],
    [0]
    ])

print(A)

sys = ct.ss(A, B, C, D)
# print(ct.poles(sys)/np.linalg.eigvals(A)(0))

Co = ct.ctrb(A,B)
nc = np.linalg.matrix_rank(Co)
print(nc)

q1 = 1
q2 = 0.1
q3 = 1
q4 = 0.1

Q = np.diag([q1, q2, q3, q4])
R = 1
K, S, E = ct.lqr(A, B, Q, R)
print(K)

sys_cl = ct.ss((A-B*K), B, C, D)
V = 0.2 * -K @ (A-B@K)
t = np.linspace(0, 10, 1000)
t, y = ct.step_response(sys, t)

# Plot the step response
plt.plot(t, y[0])
plt.xlabel("Time")
plt.ylabel("Output")
plt.title("Step Response of State-Space System")
plt.grid(True)
plt.show()