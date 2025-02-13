#!/usr/bin/env python
# coding: utf-8

# In[ ]:

# Imports
import numpy as np
import matplotlib.pyplot as plt
import integrators as intg


# Nonlinear state space form:
#  xdot = f(t, x, u), 
#   t: time
#   x: state vector
#   u: input
def f(t, x, u):
    m = 1     # in kg
    b = 0.25  # in Ns/m
    k = 1     # in N/m
    x1 = x[1]
    x2 = (-b*x[1] - k*x[0])/m
    return np.array([x1, x2])

# Time Stepping Parameters
# dt = 0.1; n = 100
# dt = 0.5; n = 100
dt = 0.1; n = 100

# Euler Integrator        
t = 0; x_euler = np.array([0, 1]); u = 0
integrator_euler = intg.Euler(dt, f)
t_history_euler = [0]
x_history_euler = [x_euler]
for i in range(n):
    
    x_euler = integrator_euler.step(t, x_euler, u)
    t = (i+1) * dt

    t_history_euler.append(t)
    x_history_euler.append(x_euler)

intg.__doc__
plt.figure()
plt.title('euler')
plt.plot(t_history_euler, x_history_euler)
plt.show()

# Heun Integrator
t = 0; x_heun = np.array([0, 1]); u = 0
integrator_heun = intg.Heun(dt, f)
t_history_heun = [0]
x_history_heun = [x_heun]
for i in range(n):
    
    x_heun = integrator_heun.step(t, x_heun, u)
    t = (i+1) * dt

    t_history_heun.append(t)
    x_history_heun.append(x_heun)

plt.figure()
plt.title('heun')
plt.plot(t_history_heun, x_history_heun)
plt.show()

# Runge Kutta 4 Integrator
t = 0; x_RK4 = np.array([0, 1]); u = 0
integrator_rk4 = intg.RK4(dt, f)
t_history_rk4 = [0]
x_history_rk4 = [x_RK4]
for i in range(n):
    
    x_RK4 = integrator_rk4.step(t, x_RK4, u)
    t = (i+1) * dt

    t_history_rk4.append(t)
    x_history_rk4.append(x_RK4)

plt.figure()
plt.title('Heun VS RK4')
plt.plot(t_history_heun, x_history_heun, label='Heun')
plt.plot(t_history_rk4, x_history_rk4, label = 'RK4')
plt.legend(["euler x", "euler x_dot", "rk4 x", "rk4 x_dot"])
plt.show()