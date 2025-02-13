import numpy as np
import matplotlib.pyplot as plt
import integrators_sonam as intg
from scipy.integrate import odeint

# Nonlinear state space form:
#  xdot = f(t, x, u), 
#   t: time
#   x: state vector
#   u: input

# Testing on Falling Ball
# def f(t, x, u):
#     return np.array([x[1], -9.81])

def f(t, x, u):
    return np.array([x[1], -x[0]-(0.25)*x[1]])

def analytical(y, t):
    x, xdot = y
    dy1dt = xdot
    dy2dt = -x-0.25*xdot
    return [dy1dt, dy2dt]
        
t = 0; x = np.array([0, 1]); u = 0
dt = 0.1; n = 50

# integrator = intg.Euler(dt, f)
integrator = intg.Heun(dt, f)
integrator1 = intg.RungeKutta4(dt, f)

t_history = [0] 
x_history = [x]
x_history1 = [x]
for i in range(n):
    
    x = integrator.step(t, x, u)
    x1 = integrator1.step(t, x, u)
    t = (i+1) * dt

    t_history.append(t)
    x_history.append(x)
    x_history1.append(x1)

t_analytic = np.linspace(0, n*dt, n)
y0 = [0, 1]
# int_analytical = odeint(analytical, y0, t_analytic)
int_analytical = (1/0.92157)*np.sin(0.992157*t_analytic)*np.exp(-0.125*t_analytic)

intg.__doc__
# plt.figure()
# plt.plot(t_analytic, int_analytical)
# plt.legend(["analytical"])
# plt.show()

xhist = np.array(x_history)

# plt.figure()
# plt.plot(t_history, x_history)
# plt.plot(t_analytic, int_analytical)
# plt.legend(["heun x", "heun x_dot", "analytical x", "analytical x_dot"])
# plt.show()

# plt.figure()
# plt.plot(t_history, xhist[:,0])
# plt.plot(t_analytic, int_analytical)
# plt.legend(["heun x","analytical x"])
# plt.show()

# plt.figure()
# plt.plot(t_history, xhist[:,1])
# plt.plot(t_analytic, int_analytical[:,1])
# plt.legend(["heun x_dot","analytical x_dot"])
# plt.show()

xhist1 = np.array(x_history1)

# plt.figure()
# plt.plot(t_history, x_history1)
# plt.plot(t_analytic, int_analytical)
# plt.legend(["rk4 x", "rk4 x_dot", "analytical x", "analytical x_dot"])
# plt.show()

plt.figure()
plt.plot(t_history, xhist1[:,0])
plt.plot(t_analytic, int_analytical)
plt.legend(["rk4 x", "analytical x"])
plt.show()

# plt.figure()
# plt.plot(t_history, xhist1[:,1])
# plt.plot(t_analytic, int_analytical[:,1])
# plt.legend(["rk4 x_dot", "analytical x_dot"])
# plt.show()