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
t = 0; x = np.array([0, 1]); u = 0
integrator_euler = intg.Euler(dt, f)
t_history_euler = [0]
x_history_euler = [x]
for i in range(n):
    
    x = integrator_euler.step(t, x, u)
    t = (i+1) * dt

    t_history_euler.append(t)
    x_history_euler.append(x)

intg.__doc__
plt.figure()
plt.title('euler')
plt.plot(t_history_euler, x_history_euler)
plt.show()

# Heun Integrator
t = 0; x = np.array([0, 1]); u = 0
integrator_heun = intg.Heun(dt, f)
t_history_heun = [0]
x_history_heun = [x]
for i in range(n):
    
    x = integrator_heun.step(t, x, u)
    t = (i+1) * dt

    t_history_heun.append(t)
    x_history_heun.append(x)

plt.figure()
plt.title('heun')
plt.plot(t_history_heun, x_history_heun)
plt.show()

# Range Kutta 4 Integrator
t = 0; x = np.array([0, 1]); u = 0
integrator_rk4 = intg.RK4(dt, f)
t_history_rk4 = [0]
x_history_rk4 = [x]
for i in range(n):
    
    x = integrator_rk4.step(t, x, u)
    t = (i+1) * dt

    t_history_rk4.append(t)
    x_history_rk4.append(x)

plt.figure()
plt.title('RK4')
plt.plot(t_history_heun, x_history_heun, label='Heun')
plt.plot(t_history_rk4, x_history_rk4, label = 'RK4')
plt.legend()
plt.show()